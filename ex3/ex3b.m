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

fig = figure(1);

subplot(2,5,[2,3,4,5])
ylim([-0.05, 0.05])
plot(time*10^6, real(on_signal), 'b')
hold on
plot(time*10^6, abs(on_signal), 'r--')
subplot(2,5,[7,8,9,10])
ylim([-0.05, 0.05])
plot(time*10^6, real(off_signal), 'b')
hold on
plot(time*10^6, abs(off_signal), 'r--')
legend('Signal', 'Envelope')

h = axes('Position',[0 0 1 1],'Visible','off');
text(.5475, .04, 'Time (\mus)')
yl = text(.2125, .43, 'Voltage V(t)');
set(yl, 'Rotation', 90)

rectangle('Position', [0, 0, 1, 1], 'EdgeColor', [0, 0, 0, 0])
rectangle('Position', [0.105, 0.6, 0.03, 0.3])
rectangle('Position', [0.075, 0.1, 0.03, 0.8])
rectangle('Position', [0.04, 0.25, 0.035, 0.03])
rectangle('Position', [0.04, 0.75, 0.035, 0.03])
annotation('doublearrow', [0.075, 0.105], [0.925, 0.925])
annotation('doublearrow', [0.105, 0.135], [0.925, 0.925])
plate_1 = text(0.085, 0.95, '10.5 mm');
plate_2 = text(0.125, 0.95, '8 mm');
set(plate_1, 'Rotation', 20)
set(plate_2, 'Rotation', 20)
text(0.04, 0.22, 'T_{OFF}')
text(0.04, 0.72, 'T_{ON}')

% Get the appropriate response - this will be the third peak when threshold
% is set to max(signal)/25.

on_threshold = max(abs(on_signal))/50;
off_threshold = max(abs(off_signal))/50;

is_response_on = logical(abs(on_signal) > on_threshold);
is_response_off = logical(abs(off_signal) > off_threshold);

this_response_on = 1;
this_response_off = 1;
for ii = 1:length(on_signal)-1
    if and(and(~is_response_on(ii), is_response_on(ii+1)), this_response_on == 5)
        % If we are at the leading edge of the third signal
        start_idx_on = ii;
    end
    if and(and(~is_response_off(ii), is_response_off(ii+1)), this_response_off == 3)
        % If we are at the leading edge of the third signal
        start_idx_off = ii;
    end
    
    if and(and(is_response_on(ii), ~is_response_on(ii+1)), this_response_on == 5)
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

len_idx = max(end_idx_on - start_idx_on, end_idx_off - start_idx_off);

% clear start_idx_on start_idx_off end_idx_on end_idx_off

fig = figure(1);

subplot(2,5,[2,3,4,5])
ylim([-0.05, 0.05])
plot(time*10^6, real(on_signal), 'b')
hold on
rectangle('Position', [time(start_idx_on)*10^6, -0.1, time(start_idx_on+len_idx)*10^6-time(start_idx_on)*10^6, 0.2], 'EdgeColor', 'none', 'FaceColor', [0.4660 0.8740 0.1880, .25])
plot(time*10^6, abs(on_signal), 'r--')
subplot(2,5,[7,8,9,10])
ylim([-0.05, 0.05])
plot(time*10^6, real(off_signal), 'b')
hold on
rectangle('Position', [time(start_idx_off)*10^6, -0.1, time(start_idx_off+len_idx)*10^6-time(start_idx_off)*10^6, 0.2], 'EdgeColor', 'none', 'FaceColor', [0.4660 0.8740 0.1880, .25])
plot(time*10^6, abs(off_signal), 'r--')
legend('Signal', 'Envelope')

h = axes('Position',[0 0 1 1],'Visible','off');
text(.5475, .04, 'Time (\mus)')
yl = text(.2125, .43, 'Voltage V(t)');
set(yl, 'Rotation', 90)

rectangle('Position', [0, 0, 1, 1], 'EdgeColor', [0, 0, 0, 0])
rectangle('Position', [0.105, 0.6, 0.03, 0.3])
rectangle('Position', [0.075, 0.1, 0.03, 0.8])
rectangle('Position', [0.04, 0.25, 0.035, 0.03])
rectangle('Position', [0.04, 0.75, 0.035, 0.03])
annotation('doublearrow', [0.075, 0.105], [0.925, 0.925])
annotation('doublearrow', [0.105, 0.135], [0.925, 0.925])
plate_1 = text(0.085, 0.95, '10.5 mm');
plate_2 = text(0.125, 0.95, '8 mm');
set(plate_1, 'Rotation', 20)
set(plate_2, 'Rotation', 20)
text(0.04, 0.22, 'T_{OFF}')
text(0.04, 0.72, 'T_{ON}')

figure(2)
plot(time(start_idx_on:start_idx_on+len_idx), real(on_signal(start_idx_on:start_idx_on+len_idx)), 'b')
hold on
plot(time(start_idx_on:start_idx_on+len_idx), abs(on_signal(start_idx_on:start_idx_on+len_idx)), 'r--')
plot(time(start_idx_off:start_idx_off+len_idx), real(off_signal(start_idx_off:start_idx_off+len_idx))+0.1, 'b')
plot(time(start_idx_off:start_idx_off+len_idx), abs(off_signal(start_idx_off:start_idx_off+len_idx))+0.1, 'r--')


% Get the frequency spectra.

refl_sig_on = on_signal(start_idx_on:start_idx_on+len_idx);
refl_sig_off = off_signal(start_idx_off:start_idx_off+len_idx);
refl_time_on = time(start_idx_on:start_idx_on+len_idx);
refl_time_off = time(start_idx_off:start_idx_off+len_idx);

fft_pts_on = 2^nextpow2(length(refl_sig_on));
fft_pts_off = 2^nextpow2(length(refl_sig_off));
refl_spec_on = fft(refl_sig_on, fft_pts_on);
refl_spec_on = refl_spec_on(1:fft_pts_on/2);
refl_spec_off = fft(refl_sig_off, fft_pts_off);
refl_spec_off = refl_spec_off(1:fft_pts_off/2);

freq_step_on = 1/(refl_time_on(end) - refl_time_on(1));
freq_on = [0 : freq_step_on : freq_step_on*(length(refl_spec_on)-1)];
freq_step_off = 1/(refl_time_off(end) - refl_time_off(1));
freq_off = [0 : freq_step_off : freq_step_off*(length(refl_spec_off)-1)];

% Calculate the theoretical reflection coefficient

z_water = rho_water * c_l_water;
z_alum = rho_alum * c_l_alum;

R_21_off = (z_alum - z_water) / (z_alum + z_water);

% Get adhesive reflection coefficient.

refl_threshold = max(abs(refl_spec_on))/100;

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

refl_spec_on_1 = refl_spec_on(start_idx:end_idx);
refl_spec_off_1 = refl_spec_off(start_idx:end_idx);
freq = freq_on(start_idx:end_idx);

R_21_on = R_21_off * real(refl_spec_on_1 ./ refl_spec_off_1);



fig = figure(3);
subplot(1,2,1)
% Only plot up to 15 MHz. Could find an automated way to do this if
% required.
plot(freq_on(1:16)*10^-6, abs(refl_spec_on(1:16)), 'b')
hold on
plot(freq_off(1:16)*10^-6, abs(refl_spec_off(1:16)), 'black')
plot([freq_on(1)*10^-6, freq_on(15)*10^-6], [refl_threshold, refl_threshold], 'r--')
box on
legend('B_1^{ON}', 'B_1^{OFF}', 'Threshold')
ylabel('Voltage |V(\omega)|')
text(.5, 2.4, '(a)')

subplot(1,2,2)
scatter(freq(2:end-1)*10^-6, R_21_on(2:end-1), 'ro')
box on
ylabel('Reflection Coefficient')
ylim([0, 1.05])
xlim([0, 20])
text(0.5, 1.0175, '(b)')

h = axes(fig,'Visible','off');
h.XLabel.Visible='on';
h.YLabel.Visible='on';
xlabel('Frequency (MHz)')