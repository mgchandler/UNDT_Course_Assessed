clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window

% Load data.

fname_on = 'joint_on_adhesive.mat';
fname_off = 'joint_off_adhesive.mat';

load(fname_on)
voltage_on = voltage;
clear time voltage

load(fname_off)
voltage_off = voltage;
clear voltage fname_on fname_off

rho_alum = 2700.0;
rho_water = 1000.0;
c_l_alum = 6320.0;
c_l_water = 1500.0;

% Filter signal.

fft_pts = 2^nextpow2(length(time));
on_spec = fft(voltage_on, fft_pts);
on_spec = on_spec(1:fft_pts/2);
off_spec = fft(voltage_off, fft_pts);
off_spec = off_spec(1:fft_pts/2);

freq_step = 1/(time(end) - time(1));
freq = [0 : freq_step : freq_step*(length(on_spec) - 1)];

window = fn_hanning(length(on_spec), 10e6/freq(end), 10e6/freq(end));
on_spec = on_spec .* window;
off_spec = off_spec .* window;

on_signal = ifft(on_spec, fft_pts);
on_signal = on_signal(1:length(voltage_on));
off_signal = ifft(off_spec, fft_pts);
off_signal = off_signal(1:length(voltage_off));

% Get the appropriate response - this will be the third peak when threshold
% is set to max(signal)/25.

on_threshold = max(abs(on_signal))/25;
off_threshold = max(abs(off_signal))/25;

is_response_on = logical(abs(on_signal) > on_threshold);
is_response_off = logical(abs(off_signal) > off_threshold);

this_response_on = 1;
this_response_off = 1;
for ii = 1:length(on_signal)-1
    if and(and(~is_response_on(ii), is_response_on(ii+1)), this_response_on == 3)
        % If we are at the leading edge of the third signal
        start_idx_on = ii;
    end
    if and(and(~is_response_off(ii), is_response_off(ii+1)), this_response_off == 3)
        % If we are at the leading edge of the third signal
        start_idx_off = ii;
    end
    
    if and(and(is_response_on(ii), ~is_response_on(ii+1)), this_response_on == 3)
       % If we are at the trailing edge of the third response
       end_idx_on = ii+1;
    end
    if and(and(is_response_off(ii), ~is_response_off(ii+1)), this_response_off == 3)
       % If we are at the trailing edge of the third response
       end_idx_off = ii+1;
    end
    
    if and(is_response_on(ii), ~is_response_on(ii+1))
       % If we are at the trailing edge of any response
       this_response_on = this_response_on + 1;
    end
    if and(is_response_off(ii), ~is_response_off(ii+1))
       % If we are at the trailing edge of any response
       this_response_off = this_response_off + 1;
    end
end

start_idx = min(start_idx_on, start_idx_off);
end_idx = max(end_idx_on, end_idx_off);

clear start_idx_on start_idx_off end_idx_on end_idx_off

% Get the frequency spectra.

refl_sig_on = on_signal(start_idx:end_idx);
refl_sig_off = off_signal(start_idx:end_idx);
refl_time = time(start_idx:end_idx);

fft_pts = 2^nextpow2(length(refl_sig_on));
refl_spec_on = fft(refl_sig_on, fft_pts);
refl_spec_on = refl_spec_on(1:fft_pts/2);
refl_spec_off = fft(refl_sig_off, fft_pts);
refl_spec_off = refl_spec_off(1:fft_pts/2);

freq_step = 1/(refl_time(end) - refl_time(1));
freq = [0 : freq_step : freq_step*(length(refl_spec_on)-1)];

% Calculate the theoretical reflection coefficient

z_water = rho_water * c_l_water;
z_alum = rho_alum * c_l_alum;

R_21_off = (z_alum - z_water) / (z_alum + z_water);

% Get adhesive reflection coefficient.

refl_threshold = max(abs(refl_spec_on))/50;

is_freq_on = logical(abs(refl_spec_on) > refl_threshold);
is_freq_off = logical(abs(refl_spec_off) > refl_threshold);

for ii = 1:length(is_freq_on)-1
    if and(~is_freq_on(ii), is_freq_on(ii+1))
        % If we are at the start of the on_freq peak
        start_idx_on = ii;
    end
    if and(~is_freq_off(ii), is_freq_off(ii+1))
        % If we are at the start of the off_freq peak
        start_idx_off = ii;
    end
    
    if and(is_freq_on(ii), ~is_freq_on(ii+1))
        % If we are at the start of the on_freq peak
        end_idx_on = ii+1;
    end
    if and(is_freq_off(ii), ~is_freq_off(ii+1))
        % If we are at the start of the off_freq peak
        end_idx_off = ii+1;
    end
end

start_idx = max(start_idx_on, start_idx_off);
end_idx = min(end_idx_on, end_idx_off);

clear start_idx_on start_idx_off end_idx_on end_idx_off

refl_spec_on = refl_spec_on(start_idx:end_idx);
refl_spec_off = refl_spec_off(start_idx:end_idx);
freq = freq(start_idx:end_idx);

R_21_on = R_21_off * real(refl_spec_on ./ refl_spec_off);

figure(1)
plot(freq, R_21_on)