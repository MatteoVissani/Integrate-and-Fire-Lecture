function [V, w ,St,Sb, I_ext] = aEIF(t,I,par)

% Author: M. Vissani, 2018

% [V w St] = aEIF(t,I,C,gL,EL,Vt,Vp,Vr,Dt,tauw,a,b,initV,initw,pflag)
%
% t     - time vector in sec. (deltaT might influence the results; smaller values earlier spikes; 0.0005)
% I     - injected current (in nA); 1. dim: time course; 2. dim: different neurons
% C     - membrane conductance (in nF)
% gL    - leak conductance (in nS)
% EL    - leak reversal (in mV), e.g., -70
% Vt    - threshold potential for spike (in mV)
% Vp    - peak potential at spike (in mV), e.g. 20
% Vr    - reset potential after spike (in mV), e.g. Vr = EL; Vr > Vt -> bursting
% Dt    - rise slope factor (in mV); sharpness of spike
% tauw  - adaptation time constant (~Ca-activated K current inactivation; in ms)
% a     - 'sub-threshold' adaptation conductance (in nS); large a = subthreshold oscillations
% b     - 'sra' current increment (in nA); large b = strong spike-frequency adaptation
% initV - initial value of membrane potential V (in mV; default = EL)
% initw - initial value of adaptation current w (in nA; default = 0)
% pflag - plotting 1 or 0
%
% V  - membrane potential (in mV)
% w  - adaptation current (in nA)
% St - spike times
%

%

%

% convert to SI units
C     = par.C*1e-9;
gL    = par.gL*1e-9;
EL    = par.EL*1e-3;
Vt    = par.Vt*1e-3;
Vp    = par.Vp*1e-3;
Vr    = par.Vr*1e-3;
Dt    = par.Dt*1e-3;
a     = par.a*1e-9;
b     = par.b*1e-9;
tauw  = par.tauw*1e-3;
initV = par.initV*1e-3;
initw = par.initw*1e-9;
abs_T = par.abs_T*1e-3;
EE = par.EE*1e-3;
tau_syn = par.tau_syn*1e-3;
init_syn = par.init_syn;
g_ext = par.g_ext*1e-9;

% time step
dt    = t(2)-t(1);    
% get number of neurons
nNeurons = par.N;


% set flag and counters to handle refractory mechanism
t_ref = zeros(1,nNeurons);
in_abs_ref= false(1,nNeurons); 

% initialize state variables
[V, w, Sb, syn,I_ext ] = deal(zeros([length(t) nNeurons]));
V(1,:)   = initV;  % membrane potential vector
w(1,:)   = initw;  % adaptation current
syn(1,:) = init_syn; % external synapse

for ii = 1 : length(t)-1
    
    % if is in refractory period --- > mantain membrane voltage and
    % increase time after last spike or reduce refractory time
    V(ii+1,in_abs_ref) = Vr;
    t_ref(in_abs_ref) = t_ref(in_abs_ref) + dt;

	% update membrane potential(only if is not refracted)
    % launch this line for no adaptation
 	%V(ii+1,~in_abs_ref) = V(ii,~in_abs_ref) + dt/C*(-gL*(V(ii,~in_abs_ref)-EL) + gL*Dt*exp((V(ii,~in_abs_ref)-Vt)/Dt)  - g_ext*syn(ii,~in_abs_ref).*(V(ii,~in_abs_ref)-EE));
    % launch this line for adaptation
    V(ii+1,~in_abs_ref) = V(ii,~in_abs_ref) + dt/C*(-gL*(V(ii,~in_abs_ref)-EL) + gL*Dt*exp((V(ii,~in_abs_ref)-Vt)/Dt) - w(ii,~in_abs_ref) - g_ext*syn(ii,~in_abs_ref).*(V(ii,~in_abs_ref)-EE));
    % launch this line for adaptation and noise
    %V(ii+1,~in_abs_ref) = V(ii,~in_abs_ref) + dt/C*(-gL*(V(ii,~in_abs_ref)-EL) + gL*Dt*exp((V(ii,~in_abs_ref)-Vt)/Dt) - w(ii,~in_abs_ref) - g_ext*syn(ii,~in_abs_ref).*(V(ii,~in_abs_ref)-EE) + sqrt(2*10^(-18))*randn());
  
    % check if the neuron will be outside its refractory time in the next
    % step
    in_abs_ref(t_ref >= abs_T) = false;
    t_ref(t_ref >= abs_T) = 0;
%     
	% check for spikes
	ixf = V(ii+1,:) > Vp;
	
	% update adaptation conductance (update always!)
	w(ii+1,:) = w(ii,:) + dt/tauw*(a*(V(ii,:)-EL) - w(ii,:)) + ixf*b;
    % update external conductance (update always!)
    syn(ii+1,:) = syn(ii,:)*(1 - dt/tau_syn) + I(ii+1,:);
	
    % set spike and reset membrane potential
	V(ii,ixf)   = Vp;
    in_abs_ref(ixf) = true;
    V(ii+1,ixf)   = Vr;
    % save binary spikes train format
	Sb(ii,ixf)  = 1; % spike bits
end

I_ext = -g_ext*syn.*(V-EE);

% get spike times
if nNeurons == 1
    St = find(Sb)*dt;
else
    St = cell([nNeurons 1]);
    for ii = 1: nNeurons
        St{ii} = find(Sb(:,ii))*dt;
    end
end

% convert back to original units
V = V*1e3;
w = w*1e9;
I_ext = I_ext*1e9;
