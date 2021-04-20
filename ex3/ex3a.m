clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window

disp('Frequency dependent attenuation of perspex'); %display the title

%INPUTS

rho_water = 1000.0;
rho_perspex = 1180.0;
rho_air = 1.2;

c_l_perspex = 2730.0;
c_l_water = 1500.0;
c_l_air = 330.0;

d = 7.8e-3;

%file name
fname = '7_8mm_thick_perspex.mat';
sample_thickness = 7.8e-3;

%load file
load(fname);

%plot the data supplied
% plot(time(5:end), voltage(5:end));
% xlabel('Time (\mus)');
% ylabel('Voltage (V)');
% text(6.75e-5, 0.2, "F")
% text(7.3e-5, 0.1, "B1")
% text(7.75e-5, 0.075, "B2")

%write your own code here to determine the frequency-dependent attenuation
%of perspex from the data supplied
%%

% Isolate signals.
%   To do this, get the envelope of the voltage signal, and get anything
%   above a defined threshold as the individual signals.

fft_pts = 2^nextpow2(length(time));
voltage_spec = fft(voltage, fft_pts);
voltage_spec = voltage_spec(1:fft_pts/2);

freq_step = 1/(time(end)-time(1));
freq = [0 : freq_step : freq_step*length(voltage_spec)-1];

voltage_signal = ifft(voltage_spec, fft_pts);
voltage_envelope = abs(voltage_signal(1:length(voltage)));

threshold = max(voltage_envelope)/100;

is_signal = logical(voltage_envelope >= threshold);

%   There are exactly three signals. Collect them separately.

individual_signals = zeros(length(voltage), 3);

start_idxs = zeros(3,1);
end_idxs = zeros(3,1);

this_signal = 1;
for ii = 1:length(is_signal)-1
    if and(~is_signal(ii), is_signal(ii+1))
        % If we are at the leading edge of a signal.
        start_idx = ii;
        start_idxs(this_signal) = start_idx;
    end
    
    if and(is_signal(ii), ~is_signal(ii+1))
        % If we are at the trailing edge of a signal.
        end_idx = ii+1;
        end_idxs(this_signal) = end_idx;
        
        % Use a square window to obtain the signal.
        individual_signals(start_idx:end_idx, this_signal) = ...
            voltage(start_idx:end_idx);
        
        this_signal = this_signal+1;
    end
end

% Calculate individual spectra
window_1 = fn_hanning(end_idxs(1)-start_idxs(1)+1, .5, .5);
spectrum_1 = fft(voltage(start_idxs(1):end_idxs(1)) .* window_1, fft_pts);
spectrum_1 = spectrum_1(1:fft_pts/2);
window_2 = fn_hanning(end_idxs(2)-start_idxs(2)+1, .5, .5);
spectrum_2 = fft(voltage(start_idxs(2):end_idxs(2)) .* window_2, fft_pts);
spectrum_2 = spectrum_2(1:fft_pts/2);
window_3 = fn_hanning(end_idxs(3)-start_idxs(3)+1, .5, .5);
spectrum_3 = fft(voltage(start_idxs(3):end_idxs(3)) .* window_3, fft_pts);
spectrum_3 = spectrum_3(1:fft_pts/2);


%plot the data supplied
figure(1)
rectangle('Position', [time(start_idxs(1)+5), -0.2, time(end_idxs(1))-time(start_idxs(1)), 0.4], 'EdgeColor', 'none', 'FaceColor', [0.4660 0.8740 0.1880, .25])
rectangle('Position', [time(start_idxs(2)+5), -0.2, time(end_idxs(2))-time(start_idxs(2)), 0.4], 'EdgeColor', 'none', 'FaceColor', [0.4660 0.8740 0.1880, .25])
rectangle('Position', [time(start_idxs(3)+5), -0.2, time(end_idxs(3))-time(start_idxs(3)), 0.4], 'EdgeColor', 'none', 'FaceColor', [0.4660 0.8740 0.1880, .25])
hold on
plot(time(5:end), real(voltage_signal(5:length(voltage))), 'b');
plot(time(5:end), voltage_envelope(5:end), 'b--')
xlabel('Time (\mus)');
ylabel('Voltage (V)');
text(6.75e-5, 0.15, "F")
text(7.3e-5, 0.075, "B1")
text(7.75e-5, 0.025, "B2")
plot([time(5), time(end)], [threshold, threshold], 'r--')
% plot([time(start_idxs(1)+5), time(start_idxs(1)+5)], [-0.2, 0.2], 'black')
% plot([time(end_idxs(1)+5), time(end_idxs(1)+5)], [-0.2, 0.2], 'black')

figure(2)
plot(freq, abs(spectrum_1))
figure(3)
plot(freq, abs(spectrum_2))
figure(4)
plot(freq, abs(spectrum_3))

freq_threshold = max(abs(spectrum_3))/5;

F_B1_spec = spectrum_2 ./ spectrum_1;
F_B1_spec(~and(abs(spectrum_1) > freq_threshold, abs(spectrum_2) > freq_threshold)) = [];
F_B1_freq = freq(and(abs(spectrum_1) > freq_threshold, abs(spectrum_2) > freq_threshold));

B1_B2_spec = spectrum_3 ./ spectrum_2;
B1_B2_spec(~and(abs(spectrum_2) > freq_threshold, abs(spectrum_3) > freq_threshold)) = [];
B1_B2_freq = freq(and(abs(spectrum_2) > freq_threshold, abs(spectrum_3) > freq_threshold));

F_B2_spec = spectrum_3 ./ spectrum_1;
F_B2_spec(~and(abs(spectrum_1) > freq_threshold, abs(spectrum_3) > freq_threshold)) = [];
F_B2_freq = freq(and(abs(spectrum_1) > freq_threshold, abs(spectrum_3) > freq_threshold));




% Calculate reflection and transmission coefficients

z_air = rho_air * c_l_air;
z_water = rho_water * c_l_water;
z_perspex = rho_perspex * c_l_perspex;

R_12 = (z_water - z_perspex) / (z_water + z_perspex);
R_21 = (z_perspex - z_water) / (z_perspex + z_water);
T_12 = 2 * z_water / (z_water + z_perspex);
T_21 = 2 * z_perspex / (z_perspex + z_water);

alpha_F_B1 = -1 / (2 * d) * log(abs(F_B1_spec * R_12 / (T_12 * R_21 * T_21)));
alpha_B1_B2 = -1 / (2 * d) * log(abs(B1_B2_spec / (R_21^2)));
alpha_F_B2 = -1 / (4 * d) * log(abs(F_B2_spec * R_12 / (T_12 * R_21^3 * T_21)));

% plot(freq, spectrum_1)
% hold on
% plot(freq, spectrum_2)

figure(5)
scatter(F_B1_freq, alpha_F_B1);
hold on
scatter(B1_B2_freq, alpha_B1_B2);
scatter(F_B2_freq, alpha_F_B2);