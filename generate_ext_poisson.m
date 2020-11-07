function ext_input=generate_ext_poisson(v_ext_eff,N,dt, steps,I_start_step, I_stop_step, I_length_step)
% ###################################
% Author:  M. Vissani, 2018
% ###################################
    ext_input = zeros(N,steps);
    ext_input(1:N,I_start_step : I_stop_step)=(rand(N,I_length_step) < ones(N,I_length_step).*(v_ext_eff*dt));
end

