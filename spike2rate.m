function rate=spike2rate(spike_count,dt,time_window,time_slide)
% ###################################
% Author:  M. Vissani, 2018
% ###################################
    slide_steps = floor(time_slide/dt);
    length_steps = floor(length(spike_count)/slide_steps);
    windows_steps= floor(floor(time_window/dt)/slide_steps);

    spike_hist_tmp = spike_count(1:(length_steps*slide_steps));
    spike_hist_tmp = reshape(spike_hist_tmp,slide_steps,length_steps);
    spike_hist = sum(spike_hist_tmp);
    kernel = zeros(1,length_steps);
    kernel(1:windows_steps) = 1;
    rate_tmp = conv(spike_hist,kernel);
    rate = rate_tmp(windows_steps:length_steps)/time_window;
end
