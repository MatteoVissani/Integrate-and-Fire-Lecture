%% Single Neuron
% ###############################
% Author: M. Vissani, 2018
% ##############################

close all
clear all
clc

h = 0.0005; % time step in s
t = 0:h:2; % time vector in s
N = 1; % number neurons layer
time_window = 50E-3;
time_slide = 5E-3;

% Stimulus Parameter

% Time Parameter
I_start = 0.5; % stimulus onset in s
I_start_step = floor(I_start/h);
I_stop = 1; % stimulus offset in s
I_stop_step = floor(I_stop/h);
I_length = I_stop - I_start; % stimulus length
I_length_step = I_stop_step - I_start_step + 1;


% Frequency Parameter
v_ext = 300*ones(1,N)'; % stimulus frequency

steps = numel(t); % number of steps
I = generate_ext_poisson(v_ext,N,h, steps, I_start_step, I_stop_step, I_length_step);

% spike2rate convert the spikes 1/0 binary format in firing rate Hz
rate_ext = spike2rate(I,h,time_window,time_slide); % plot average input stimulus
figure
x=0:time_slide:2;
x=x(1:length(rate_ext));
plot(x,rate_ext,'r');
xlabel('time(s)');
ylabel('firing rate(Hz)');
title(' average input firing rate')


% Network and Neuron Parameters
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
% abs_T -  absolute refractory period (in ms; default = 2) 
%
% V  - membrane potential (in mV)
% w  - adaptation current (in nA)
% St - spike times
%
par = struct();
par.C = 0.45;
par.gL = 25;
par.EL = -70.6;
par.Vt = -50.4;% -50.4
par.Vp = 20;%20
par.Vr = -47.4;
par.Dt = 2;
par.tauw = 144;%144
par.a = 4;%4
par.b = 0.0805;%0.0805
par.initV = -70.6;
par.initw = 0;
par.tau_syn = 10;
par.init_syn = 0;
par.g_ext = 35; %35
par.EE = 0;%0
par.N = N;
par.abs_T = 2; 
I = I';

% Build Neurons and Dynamics ( see the function for parameters)

[V, w, St,Sb, I_ext] = aEIF(t,I,par);
rate = spike2rate(Sb,h,time_window,time_slide);

% % plot
figure, set(gcf,'Color',[1 1 1])
subplot(6,1,1:3)
plot(t,V,'-k')
title('AdEx LIF neuron')
ylabel('V (mV)')
axis([t([1 end]) -100 30])
subplot(6,1,4:5)
plot(t,w,'-k')
ylabel('w (nA)')
axis([t([1 end]) min(w)-0.02 max(w)+0.02])
subplot(6,1,6)
plot(t,I_ext,'-k')
xlabel('t (s)'); ylabel('Input (nA)')
axis([t([1 end]) min(I_ext)-0.02 max(I_ext)+0.02])

figure

plot(rate)
xlabel('t (s)'); ylabel('Firing Rate [Hz] ')
title(' Neuron activity ')

%% animated plot
% f = figure('renderer','painters');
% anim_V = animatedline('Color','k','LineWidth',1);
% set(gca,'Xlim',[t(1) t(end)],'Ylim',[-100 30])
% xlabel('Time [s]')
% ylabel('Voltage [mV]')
% title('In-silico Neuron')
% grid on
% grid minor
% gif_name = 'Insilico_Neuron.gif';
% dsample_factor = 10; % to improve speed. Only for graphical purposes
% 
% for ii = 1:steps
%     addpoints(anim_V,t(ii),V(ii));
%     drawnow
%     
%     if mod(ii,dsample_factor) == 0 % downsampling condition
%        % Capture the plot as an image 
%        frame = getframe(f); 
%        im = frame2im(frame); 
%        [imind,cm] = rgb2ind(im,256); 
%         % Write to the GIF File 
%         if ii == 1 
%            imwrite(imind,cm,gif_name,'gif', 'Loopcount',inf); 
%         else 
%            imwrite(imind,cm,gif_name,'gif','WriteMode','append','DelayTime',0.1); 
%         end 
%     end
% end