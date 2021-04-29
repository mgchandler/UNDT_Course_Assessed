clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window

%% This is a new code which only focusses on TFM to identify optimal imaging
% parameters for the calibration sample.

%% Inputs

num_els = 64;
centre_freq_1 = 5e6;
el_sep = .05e-3;
v_L = 1500.0;
el_pitch_1 = 0.15e-3;
backwall_dist = 20.0e-3;
grid_pts = 351;
scat_coords = [[-2.25e-3,  5.00e-3]; ...
               [-1.00e-3,  7.50e-3]; ...
               [ 0.25e-3, 10.00e-3]; ...
               [ 1.50e-3, 12.50e-3]; ...
               [ 2.75e-3, 15.00e-3]; ...
               [ 0.00e-3,  5.00e-3]; ...
               [ 0.00e-3,  7.50e-3]; ...
               [ 0.00e-3, 10.00e-3]; ...
               [ 0.00e-3, 12.50e-3]; ...
               [ 0.00e-3, 15.00e-3]];
           
% Miscellaneous extra params.
num_cycles = 5;
oversample = 10; % How often per cycle the signal is sampled.

%% Parameters

lambda = v_L / centre_freq_1;

probe_coords_1 = zeros(num_els, 2);
probe_coords_1(:, 1) = linspace(0, el_pitch_1*(num_els-1), num_els);
probe_coords_1(:, 1) = probe_coords_1(:, 1) - mean(probe_coords_1(:, 1));

[tx, rx] = meshgrid([1:num_els]);
tx = reshape(tx, length(tx)^2, 1);
rx = reshape(rx, length(rx)^2, 1);

%% Get input time signal and freq spectrum

max_time = 3 * sqrt(2 * backwall_dist^2) / v_L;
dt = 1 / centre_freq_1 / oversample;
time_1 = [0 : dt : max_time];

half_pulse = num_cycles / (2 * centre_freq_1);

input_signal = sin(2 * pi * centre_freq_1 * time_1)' .* ...
               fn_hanning(length(time_1), half_pulse / max_time, half_pulse / max_time);
           
% Bring peak of input pulse to time = 0.
time_1 = time_1 - half_pulse;
           
           
           
fft_pts = 2^nextpow2(length(time_1));

df = 1 / max_time;
freq = [0 : df : df*(fft_pts/2-1)];
omega = 2 * pi * freq;

input_spectrum = fft(input_signal, fft_pts);
input_spectrum = input_spectrum(1:fft_pts/2);

%% Do ray tracing

scat_dists_tx = sqrt((probe_coords_1(tx, 1) - scat_coords(:, 1)').^2 + ...
                      (probe_coords_1(tx, 2) - scat_coords(:, 2)').^2);
scat_dists_rx = sqrt((probe_coords_1(rx, 1) - scat_coords(:, 1)').^2 + ...
                     (probe_coords_1(rx, 2) - scat_coords(:, 2)').^2);
scat_times = (scat_dists_tx + scat_dists_rx) / v_L;
                      
bw_dists = (sqrt((probe_coords_1(tx, 1)' - probe_coords_1(rx, 1)').^2 + ...
                 (2*(backwall_dist - probe_coords_1(tx, 2)')).^2));
bw_times = bw_dists / v_L;

% Get angles for directivity.
scat_theta_tx = acos((scat_coords(:, 2)' - probe_coords_1(tx, 2)) ./ ...
                     sqrt((probe_coords_1(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_1(tx, 2) - scat_coords(:, 2)').^2));
scat_theta_rx = acos((scat_coords(:, 2)' - probe_coords_1(rx, 2)) ./ ...
                     sqrt((probe_coords_1(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_1(tx, 2) - scat_coords(:, 2)').^2));
                      
bw_theta_tx = acos((2*(backwall_dist - probe_coords_1(tx, 2)')) ./ bw_dists);
bw_theta_rx = acos((2*(backwall_dist - probe_coords_1(rx, 2)')) ./ bw_dists);

%% Propagate signal

% Get scatterer and backwall signal amplitudes.

% Directivity includes both transmit and receive paths. N.B. this function
% takes into account that in Matlab, sinc(x) := sin(πx)/(πx)
scat_dir = (el_pitch_1 - el_sep)^2 * ...
           sinc((el_pitch_1 - el_sep) / lambda * sin(scat_theta_tx)) .* ...
           sinc((el_pitch_1 - el_sep) / lambda * sin(scat_theta_rx));
% Amplitude includes scattering, directivity and beam spreading.
scat_amp = 0.01 * scat_dir ./ (sqrt(scat_dists_tx) .* sqrt(scat_dists_rx));
                               
bw_dir = (el_pitch_1 - el_sep)^2 * ...
         sinc((el_pitch_1 - el_sep) / lambda * sin(bw_theta_tx)) .* ...
         sinc((el_pitch_1 - el_sep) / lambda * sin(bw_theta_rx));
bw_amp = bw_dir ./ sqrt(bw_dists); % Perfect reflector => R = 1.



output_spectra = 0;
% Scatterers
for ss = 1:size(scat_coords, 1)
    output_spectra = output_spectra + ( ...
        spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
        exp(-1i * omega' * scat_times(:, ss)') * ...
        spdiags(scat_amp(:, ss), 0, length(tx), length(tx)) ...
    );
end

% Backwall
output_spectra = output_spectra + ( ...
    spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
    exp(-1i * omega' * bw_times(:)') * ...
    spdiags(bw_amp(:), 0, length(tx), length(tx)) ...
);

FMC_data_1 = ifft(output_spectra, fft_pts, 1);
FMC_data_1 = FMC_data_1(1:length(time_1), :); 
FMC_dt = 1 / (fft_pts * abs(freq(2) - freq(1)));
FMC_time_1 = [0 : FMC_dt : FMC_dt * (length(time_1)-1)];

% End of simulation of wave.

%% Imaging - TFM

x = linspace(- backwall_dist / 2 - 1e-3, backwall_dist / 2 + 1e-3, grid_pts);
z = linspace(-1e-3, backwall_dist + 1e-3, grid_pts);
% x = linspace(-3.5e-3, 3.5e-3, grid_pts);
% z = linspace(7e-3, 14e-3, grid_pts);
[X, Z] = meshgrid(x, z);

Focus_lookup_times = (sqrt((probe_coords_1(tx, 1) - X(:)').^2 + ...
                           (probe_coords_1(tx, 2) - Z(:)').^2) + ...
                      sqrt((probe_coords_1(rx, 1) - X(:)').^2 + ...
                           (probe_coords_1(rx, 2) - Z(:)').^2) ...
    ) / v_L;

% TFM lookup times are identical to focussed lookup times. Difference is
% that we focus everywhere, instead of just below the aperture.

Im_TFM = 0;
tau_TFM = reshape(Focus_lookup_times, num_els^2, size(X, 1), size(X, 2));
for tr = 1:num_els^2
    Im_TFM = Im_TFM + ...
        interp1(FMC_time_1, FMC_data_1(:, tr), squeeze(tau_TFM(tr, :, :)), 'linear', 0);
end

max_ = max(abs(Im_TFM), [], 'all');
dB_TFM_1 = 20 * log10(abs(Im_TFM) ./ max_);
disp(el_pitch_1)



%% Inputs

num_els = 64;
centre_freq_2 = 5.25e6;
el_sep = .05e-3;
v_L = 1500.0;
el_pitch_2 = 0.15e-3;
backwall_dist = 20.0e-3;
grid_pts = 351;
scat_coords = [[-2.25e-3,  5.00e-3]; ...
               [-1.00e-3,  7.50e-3]; ...
               [ 0.25e-3, 10.00e-3]; ...
               [ 1.50e-3, 12.50e-3]; ...
               [ 2.75e-3, 15.00e-3]; ...
               [ 0.00e-3,  5.00e-3]; ...
               [ 0.00e-3,  7.50e-3]; ...
               [ 0.00e-3, 10.00e-3]; ...
               [ 0.00e-3, 12.50e-3]; ...
               [ 0.00e-3, 15.00e-3]];
           
% Miscellaneous extra params.
num_cycles = 5;
oversample = 10; % How often per cycle the signal is sampled.

%% Parameters

lambda = v_L / centre_freq_2;

probe_coords_2 = zeros(num_els, 2);
probe_coords_2(:, 1) = linspace(0, el_pitch_2*(num_els-1), num_els);
probe_coords_2(:, 1) = probe_coords_2(:, 1) - mean(probe_coords_2(:, 1));

[tx, rx] = meshgrid([1:num_els]);
tx = reshape(tx, length(tx)^2, 1);
rx = reshape(rx, length(rx)^2, 1);

%% Get input time signal and freq spectrum

max_time = 3 * sqrt(2 * backwall_dist^2) / v_L;
dt = 1 / centre_freq_2 / oversample;
time_2 = [0 : dt : max_time];

half_pulse = num_cycles / (2 * centre_freq_2);

input_signal = sin(2 * pi * centre_freq_2 * time_2)' .* ...
               fn_hanning(length(time_2), half_pulse / max_time, half_pulse / max_time);
           
% Bring peak of input pulse to time = 0.
time_2 = time_2 - half_pulse;
           
           
           
fft_pts = 2^nextpow2(length(time_2));

df = 1 / max_time;
freq = [0 : df : df*(fft_pts/2-1)];
omega = 2 * pi * freq;

input_spectrum = fft(input_signal, fft_pts);
input_spectrum = input_spectrum(1:fft_pts/2);

%% Do ray tracing

scat_dists_tx = sqrt((probe_coords_2(tx, 1) - scat_coords(:, 1)').^2 + ...
                      (probe_coords_2(tx, 2) - scat_coords(:, 2)').^2);
scat_dists_rx = sqrt((probe_coords_2(rx, 1) - scat_coords(:, 1)').^2 + ...
                     (probe_coords_2(rx, 2) - scat_coords(:, 2)').^2);
scat_times = (scat_dists_tx + scat_dists_rx) / v_L;
                      
bw_dists = (sqrt((probe_coords_2(tx, 1)' - probe_coords_2(rx, 1)').^2 + ...
                 (2*(backwall_dist - probe_coords_2(tx, 2)')).^2));
bw_times = bw_dists / v_L;

% Get angles for directivity.
scat_theta_tx = acos((scat_coords(:, 2)' - probe_coords_2(tx, 2)) ./ ...
                     sqrt((probe_coords_2(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_2(tx, 2) - scat_coords(:, 2)').^2));
scat_theta_rx = acos((scat_coords(:, 2)' - probe_coords_2(rx, 2)) ./ ...
                     sqrt((probe_coords_2(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_2(tx, 2) - scat_coords(:, 2)').^2));
                      
bw_theta_tx = acos((2*(backwall_dist - probe_coords_2(tx, 2)')) ./ bw_dists);
bw_theta_rx = acos((2*(backwall_dist - probe_coords_2(rx, 2)')) ./ bw_dists);

%% Propagate signal

% Get scatterer and backwall signal amplitudes.

% Directivity includes both transmit and receive paths. N.B. this function
% takes into account that in Matlab, sinc(x) := sin(πx)/(πx)
scat_dir = (el_pitch_2 - el_sep)^2 * ...
           sinc((el_pitch_2 - el_sep) / lambda * sin(scat_theta_tx)) .* ...
           sinc((el_pitch_2 - el_sep) / lambda * sin(scat_theta_rx));
% Amplitude includes scattering, directivity and beam spreading.
scat_amp = 0.01 * scat_dir ./ (sqrt(scat_dists_tx) .* sqrt(scat_dists_rx));
                               
bw_dir = (el_pitch_2 - el_sep)^2 * ...
         sinc((el_pitch_2 - el_sep) / lambda * sin(bw_theta_tx)) .* ...
         sinc((el_pitch_2 - el_sep) / lambda * sin(bw_theta_rx));
bw_amp = bw_dir ./ sqrt(bw_dists); % Perfect reflector => R = 1.



output_spectra = 0;
% Scatterers
for ss = 1:size(scat_coords, 1)
    output_spectra = output_spectra + ( ...
        spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
        exp(-1i * omega' * scat_times(:, ss)') * ...
        spdiags(scat_amp(:, ss), 0, length(tx), length(tx)) ...
    );
end

% Backwall
output_spectra = output_spectra + ( ...
    spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
    exp(-1i * omega' * bw_times(:)') * ...
    spdiags(bw_amp(:), 0, length(tx), length(tx)) ...
);

FMC_data_2 = ifft(output_spectra, fft_pts, 1);
FMC_data_2 = FMC_data_2(1:length(time_2), :); 
FMC_dt = 1 / (fft_pts * abs(freq(2) - freq(1)));
FMC_time_2 = [0 : FMC_dt : FMC_dt * (length(time_2)-1)];

% End of simulation of wave.

%% Imaging - TFM

x = linspace(- backwall_dist / 2 - 1e-3, backwall_dist / 2 + 1e-3, grid_pts);
z = linspace(-1e-3, backwall_dist + 1e-3, grid_pts);
% x = linspace(-3.5e-3, 3.5e-3, grid_pts);
% z = linspace(7e-3, 14e-3, grid_pts);
[X, Z] = meshgrid(x, z);

Focus_lookup_times = (sqrt((probe_coords_2(tx, 1) - X(:)').^2 + ...
                           (probe_coords_2(tx, 2) - Z(:)').^2) + ...
                      sqrt((probe_coords_2(rx, 1) - X(:)').^2 + ...
                           (probe_coords_2(rx, 2) - Z(:)').^2) ...
    ) / v_L;

% TFM lookup times are identical to focussed lookup times. Difference is
% that we focus everywhere, instead of just below the aperture.

Im_TFM = 0;
tau_TFM = reshape(Focus_lookup_times, num_els^2, size(X, 1), size(X, 2));
for tr = 1:num_els^2
    Im_TFM = Im_TFM + ...
        interp1(FMC_time_2, FMC_data_2(:, tr), squeeze(tau_TFM(tr, :, :)), 'linear', 0);
end

max_ = max(abs(Im_TFM), [], 'all');
dB_TFM_2 = 20 * log10(abs(Im_TFM) ./ max_);
disp(el_pitch_2)



%% Inputs

num_els = 64;
centre_freq_3 = 5.5e6;
el_sep = .05e-3;
v_L = 1500.0;
el_pitch_3 = 0.15e-3;
backwall_dist = 20.0e-3;
grid_pts = 351;
scat_coords = [[-2.25e-3,  5.00e-3]; ...
               [-1.00e-3,  7.50e-3]; ...
               [ 0.25e-3, 10.00e-3]; ...
               [ 1.50e-3, 12.50e-3]; ...
               [ 2.75e-3, 15.00e-3]; ...
               [ 0.00e-3,  5.00e-3]; ...
               [ 0.00e-3,  7.50e-3]; ...
               [ 0.00e-3, 10.00e-3]; ...
               [ 0.00e-3, 12.50e-3]; ...
               [ 0.00e-3, 15.00e-3]];
           
% Miscellaneous extra params.
num_cycles = 5;
oversample = 10; % How often per cycle the signal is sampled.

%% Parameters

lambda = v_L / centre_freq_3;

probe_coords_3 = zeros(num_els, 2);
probe_coords_3(:, 1) = linspace(0, el_pitch_3*(num_els-1), num_els);
probe_coords_3(:, 1) = probe_coords_3(:, 1) - mean(probe_coords_3(:, 1));

[tx, rx] = meshgrid([1:num_els]);
tx = reshape(tx, length(tx)^2, 1);
rx = reshape(rx, length(rx)^2, 1);

%% Get input time signal and freq spectrum

max_time = 3 * sqrt(2 * backwall_dist^2) / v_L;
dt = 1 / centre_freq_3 / oversample;
time_3 = [0 : dt : max_time];

half_pulse = num_cycles / (2 * centre_freq_3);

input_signal = sin(2 * pi * centre_freq_3 * time_3)' .* ...
               fn_hanning(length(time_3), half_pulse / max_time, half_pulse / max_time);
           
% Bring peak of input pulse to time = 0.
time_3 = time_3 - half_pulse;
           
           
           
fft_pts = 2^nextpow2(length(time_3));

df = 1 / max_time;
freq = [0 : df : df*(fft_pts/2-1)];
omega = 2 * pi * freq;

input_spectrum = fft(input_signal, fft_pts);
input_spectrum = input_spectrum(1:fft_pts/2);

%% Do ray tracing

scat_dists_tx = sqrt((probe_coords_3(tx, 1) - scat_coords(:, 1)').^2 + ...
                      (probe_coords_3(tx, 2) - scat_coords(:, 2)').^2);
scat_dists_rx = sqrt((probe_coords_3(rx, 1) - scat_coords(:, 1)').^2 + ...
                     (probe_coords_3(rx, 2) - scat_coords(:, 2)').^2);
scat_times = (scat_dists_tx + scat_dists_rx) / v_L;
                      
bw_dists = (sqrt((probe_coords_3(tx, 1)' - probe_coords_3(rx, 1)').^2 + ...
                 (2*(backwall_dist - probe_coords_3(tx, 2)')).^2));
bw_times = bw_dists / v_L;

% Get angles for directivity.
scat_theta_tx = acos((scat_coords(:, 2)' - probe_coords_3(tx, 2)) ./ ...
                     sqrt((probe_coords_3(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_3(tx, 2) - scat_coords(:, 2)').^2));
scat_theta_rx = acos((scat_coords(:, 2)' - probe_coords_3(rx, 2)) ./ ...
                     sqrt((probe_coords_3(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_3(tx, 2) - scat_coords(:, 2)').^2));
                      
bw_theta_tx = acos((2*(backwall_dist - probe_coords_3(tx, 2)')) ./ bw_dists);
bw_theta_rx = acos((2*(backwall_dist - probe_coords_3(rx, 2)')) ./ bw_dists);

%% Propagate signal

% Get scatterer and backwall signal amplitudes.

% Directivity includes both transmit and receive paths. N.B. this function
% takes into account that in Matlab, sinc(x) := sin(πx)/(πx)
scat_dir = (el_pitch_3 - el_sep)^2 * ...
           sinc((el_pitch_3 - el_sep) / lambda * sin(scat_theta_tx)) .* ...
           sinc((el_pitch_3 - el_sep) / lambda * sin(scat_theta_rx));
% Amplitude includes scattering, directivity and beam spreading.
scat_amp = 0.01 * scat_dir ./ (sqrt(scat_dists_tx) .* sqrt(scat_dists_rx));
                               
bw_dir = (el_pitch_3 - el_sep)^2 * ...
         sinc((el_pitch_3 - el_sep) / lambda * sin(bw_theta_tx)) .* ...
         sinc((el_pitch_3 - el_sep) / lambda * sin(bw_theta_rx));
bw_amp = bw_dir ./ sqrt(bw_dists); % Perfect reflector => R = 1.



output_spectra = 0;
% Scatterers
for ss = 1:size(scat_coords, 1)
    output_spectra = output_spectra + ( ...
        spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
        exp(-1i * omega' * scat_times(:, ss)') * ...
        spdiags(scat_amp(:, ss), 0, length(tx), length(tx)) ...
    );
end

% Backwall
output_spectra = output_spectra + ( ...
    spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
    exp(-1i * omega' * bw_times(:)') * ...
    spdiags(bw_amp(:), 0, length(tx), length(tx)) ...
);

FMC_data_3 = ifft(output_spectra, fft_pts, 1);
FMC_data_3 = FMC_data_3(1:length(time_3), :); 
FMC_dt = 1 / (fft_pts * abs(freq(2) - freq(1)));
FMC_time_3 = [0 : FMC_dt : FMC_dt * (length(time_3)-1)];

% End of simulation of wave.

%% Imaging - TFM

x = linspace(- backwall_dist / 2 - 1e-3, backwall_dist / 2 + 1e-3, grid_pts);
z = linspace(-1e-3, backwall_dist + 1e-3, grid_pts);
% x = linspace(-3.5e-3, 3.5e-3, grid_pts);
% z = linspace(7e-3, 14e-3, grid_pts);
[X, Z] = meshgrid(x, z);

Focus_lookup_times = (sqrt((probe_coords_3(tx, 1) - X(:)').^2 + ...
                           (probe_coords_3(tx, 2) - Z(:)').^2) + ...
                      sqrt((probe_coords_3(rx, 1) - X(:)').^2 + ...
                           (probe_coords_3(rx, 2) - Z(:)').^2) ...
    ) / v_L;

% TFM lookup times are identical to focussed lookup times. Difference is
% that we focus everywhere, instead of just below the aperture.

Im_TFM = 0;
tau_TFM = reshape(Focus_lookup_times, num_els^2, size(X, 1), size(X, 2));
for tr = 1:num_els^2
    Im_TFM = Im_TFM + ...
        interp1(FMC_time_3, FMC_data_3(:, tr), squeeze(tau_TFM(tr, :, :)), 'linear', 0);
end

max_ = max(abs(Im_TFM), [], 'all');
dB_TFM_3 = 20 * log10(abs(Im_TFM) ./ max_);
disp(el_pitch_3)



%% Inputs

num_els = 64;
centre_freq_4 = 5.75e6;
el_sep = .05e-3;
v_L = 1500.0;
el_pitch_4 = 0.15e-3;
backwall_dist = 20.0e-3;
grid_pts = 351;
scat_coords = [[-2.25e-3,  5.00e-3]; ...
               [-1.00e-3,  7.50e-3]; ...
               [ 0.25e-3, 10.00e-3]; ...
               [ 1.50e-3, 12.50e-3]; ...
               [ 2.75e-3, 15.00e-3]; ...
               [ 0.00e-3,  5.00e-3]; ...
               [ 0.00e-3,  7.50e-3]; ...
               [ 0.00e-3, 10.00e-3]; ...
               [ 0.00e-3, 12.50e-3]; ...
               [ 0.00e-3, 15.00e-3]];
           
% Miscellaneous extra params.
num_cycles = 5;
oversample = 10; % How often per cycle the signal is sampled.

%% Parameters

lambda = v_L / centre_freq_4;

probe_coords_4 = zeros(num_els, 2);
probe_coords_4(:, 1) = linspace(0, el_pitch_4*(num_els-1), num_els);
probe_coords_4(:, 1) = probe_coords_4(:, 1) - mean(probe_coords_4(:, 1));

[tx, rx] = meshgrid([1:num_els]);
tx = reshape(tx, length(tx)^2, 1);
rx = reshape(rx, length(rx)^2, 1);

%% Get input time signal and freq spectrum

max_time = 3 * sqrt(2 * backwall_dist^2) / v_L;
dt = 1 / centre_freq_4 / oversample;
time_4 = [0 : dt : max_time];

half_pulse = num_cycles / (2 * centre_freq_4);

input_signal = sin(2 * pi * centre_freq_4 * time_4)' .* ...
               fn_hanning(length(time_4), half_pulse / max_time, half_pulse / max_time);
           
% Bring peak of input pulse to time = 0.
time_4 = time_4 - half_pulse;
           
           
           
fft_pts = 2^nextpow2(length(time_4));

df = 1 / max_time;
freq = [0 : df : df*(fft_pts/2-1)];
omega = 2 * pi * freq;

input_spectrum = fft(input_signal, fft_pts);
input_spectrum = input_spectrum(1:fft_pts/2);

%% Do ray tracing

scat_dists_tx = sqrt((probe_coords_4(tx, 1) - scat_coords(:, 1)').^2 + ...
                      (probe_coords_4(tx, 2) - scat_coords(:, 2)').^2);
scat_dists_rx = sqrt((probe_coords_4(rx, 1) - scat_coords(:, 1)').^2 + ...
                     (probe_coords_4(rx, 2) - scat_coords(:, 2)').^2);
scat_times = (scat_dists_tx + scat_dists_rx) / v_L;
                      
bw_dists = (sqrt((probe_coords_4(tx, 1)' - probe_coords_4(rx, 1)').^2 + ...
                 (2*(backwall_dist - probe_coords_4(tx, 2)')).^2));
bw_times = bw_dists / v_L;

% Get angles for directivity.
scat_theta_tx = acos((scat_coords(:, 2)' - probe_coords_4(tx, 2)) ./ ...
                     sqrt((probe_coords_4(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_4(tx, 2) - scat_coords(:, 2)').^2));
scat_theta_rx = acos((scat_coords(:, 2)' - probe_coords_4(rx, 2)) ./ ...
                     sqrt((probe_coords_4(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords_4(tx, 2) - scat_coords(:, 2)').^2));
                      
bw_theta_tx = acos((2*(backwall_dist - probe_coords_4(tx, 2)')) ./ bw_dists);
bw_theta_rx = acos((2*(backwall_dist - probe_coords_4(rx, 2)')) ./ bw_dists);

%% Propagate signal

% Get scatterer and backwall signal amplitudes.

% Directivity includes both transmit and receive paths. N.B. this function
% takes into account that in Matlab, sinc(x) := sin(πx)/(πx)
scat_dir = (el_pitch_4 - el_sep)^2 * ...
           sinc((el_pitch_4 - el_sep) / lambda * sin(scat_theta_tx)) .* ...
           sinc((el_pitch_4 - el_sep) / lambda * sin(scat_theta_rx));
% Amplitude includes scattering, directivity and beam spreading.
scat_amp = 0.01 * scat_dir ./ (sqrt(scat_dists_tx) .* sqrt(scat_dists_rx));
                               
bw_dir = (el_pitch_4 - el_sep)^2 * ...
         sinc((el_pitch_4 - el_sep) / lambda * sin(bw_theta_tx)) .* ...
         sinc((el_pitch_4 - el_sep) / lambda * sin(bw_theta_rx));
bw_amp = bw_dir ./ sqrt(bw_dists); % Perfect reflector => R = 1.



output_spectra = 0;
% Scatterers
for ss = 1:size(scat_coords, 1)
    output_spectra = output_spectra + ( ...
        spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
        exp(-1i * omega' * scat_times(:, ss)') * ...
        spdiags(scat_amp(:, ss), 0, length(tx), length(tx)) ...
    );
end

% Backwall
output_spectra = output_spectra + ( ...
    spdiags(input_spectrum, 0, length(freq), length(freq)) * ...
    exp(-1i * omega' * bw_times(:)') * ...
    spdiags(bw_amp(:), 0, length(tx), length(tx)) ...
);

FMC_data_4 = ifft(output_spectra, fft_pts, 1);
FMC_data_4 = FMC_data_4(1:length(time_4), :); 
FMC_dt = 1 / (fft_pts * abs(freq(2) - freq(1)));
FMC_time_4 = [0 : FMC_dt : FMC_dt * (length(time_4)-1)];

% End of simulation of wave.

%% Imaging - TFM

x = linspace(- backwall_dist / 2 - 1e-3, backwall_dist / 2 + 1e-3, grid_pts);
z = linspace(-1e-3, backwall_dist + 1e-3, grid_pts);
% x = linspace(-3.5e-3, 3.5e-3, grid_pts);
% z = linspace(7e-3, 14e-3, grid_pts);
[X, Z] = meshgrid(x, z);

Focus_lookup_times = (sqrt((probe_coords_4(tx, 1) - X(:)').^2 + ...
                           (probe_coords_4(tx, 2) - Z(:)').^2) + ...
                      sqrt((probe_coords_4(rx, 1) - X(:)').^2 + ...
                           (probe_coords_4(rx, 2) - Z(:)').^2) ...
    ) / v_L;

% TFM lookup times are identical to focussed lookup times. Difference is
% that we focus everywhere, instead of just below the aperture.

Im_TFM = 0;
tau_TFM = reshape(Focus_lookup_times, num_els^2, size(X, 1), size(X, 2));
for tr = 1:num_els^2
    Im_TFM = Im_TFM + ...
        interp1(FMC_time_4, FMC_data_4(:, tr), squeeze(tau_TFM(tr, :, :)), 'linear', 0);
end

max_ = max(abs(Im_TFM), [], 'all');
dB_TFM_4 = 20 * log10(abs(Im_TFM) ./ max_);
disp(el_pitch_4)



%% Plotting

% figure(1)
% t1 = tiledlayout(2, 3);
% nexttile;
% imagesc(FMC_time_1*10^6, [1:num_els^2], abs(FMC_data_1'))
% box on
% text(1, 120, '(a)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_1*10^-6, el_pitch_1*10^3))
% 
% nexttile;
% imagesc(FMC_time_2*10^6, [1:num_els^2], abs(FMC_data_2'))
% box on
% text(1, 120, '(b)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_2*10^-6, el_pitch_2*10^3))
% 
% nexttile;
% imagesc(FMC_time_3*10^6, [1:num_els^2], abs(FMC_data_3'))
% box on
% text(1, 120, '(c)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_3*10^-6, el_pitch_3*10^3))
% 
% nexttile;
% imagesc(FMC_time_4*10^6, [1:num_els^2], abs(FMC_data_4'))
% box on
% text(1, 120, '(d)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_4*10^-6, el_pitch_4*10^3))
% 
% nexttile;
% imagesc(FMC_time_5*10^6, [1:num_els^2], abs(FMC_data_5'))
% box on
% text(1, 120, '(e)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_5*10^-6, el_pitch_5*10^3))
% 
% nexttile;
% imagesc(FMC_time_6*10^6, [1:num_els^2], abs(FMC_data_6'))
% box on
% text(1, 120, '(f)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_6*10^-6, el_pitch_6*10^3))
% 
% xlabel(t1, 'Time (μs)')
% ylabel(t1, 'Tx-Rx pair')
% 
% cb = colorbar();
% cb.Label.String = 'dB';
% cb.Layout.Tile = 'east';

figure(1)
t2 = tiledlayout(2, 4);
nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_1, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_1(:,1)*10^3, probe_coords_1(:,2)*10^3, 'wo')
rectangle('Position', [-1.5, 9, 3, 3], 'EdgeColor', 'r')
% xlim([-1.5, 1.5])
% ylim([9, 12])
text(-10, 19, '(a)', 'Color', 'w')
title(sprintf('f = %3.2fMHz', centre_freq_1*10^-6))

nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_2, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_2(:,1)*10^3, probe_coords_2(:,2)*10^3, 'wo')
rectangle('Position', [-1.5, 9, 3, 3], 'EdgeColor', 'r')
% xlim([-1.5, 1.5])
% ylim([9, 12])
text(-10, 19, '(b)', 'Color', 'w')
title(sprintf('f = %3.2fMHz', centre_freq_2*10^-6))

nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_3, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_3(:,1)*10^3, probe_coords_3(:,2)*10^3, 'wo')
rectangle('Position', [-1.5, 9, 3, 3], 'EdgeColor', 'r')
% xlim([-1.5, 1.5])
% ylim([9, 12])
text(-10, 19, '(c)', 'Color', 'w')
title(sprintf('f = %3.2fMHz', centre_freq_3*10^-6))

nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_4, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_4(:,1)*10^3, probe_coords_4(:,2)*10^3, 'wo')
rectangle('Position', [-1.5, 9, 3, 3], 'EdgeColor', 'r')
% xlim([-1.5, 1.5])
% ylim([9, 12])
text(-10, 19, '(d)', 'Color', 'w')
title(sprintf('f = %3.2fMHz', centre_freq_4*10^-6))

dB_TFM_1(dB_TFM_1 < -10) = -40;
dB_TFM_2(dB_TFM_2 < -10) = -40;
dB_TFM_3(dB_TFM_3 < -10) = -40;
dB_TFM_4(dB_TFM_4 < -10) = -40;

nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_1, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_1(:,1)*10^3, probe_coords_1(:,2)*10^3, 'wo')
xlim([-1.5, 1.5])
ylim([9, 12])
text(-1.4, 11.7, '(e)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_1*10^-6, el_pitch_1*10^3))

nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_2, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_2(:,1)*10^3, probe_coords_2(:,2)*10^3, 'wo')
xlim([-1.5, 1.5])
ylim([9, 12])
text(-1.4, 11.7, '(f)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_2*10^-6, el_pitch_2*10^3))

nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_3, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_3(:,1)*10^3, probe_coords_3(:,2)*10^3, 'wo')
xlim([-1.5, 1.5])
ylim([9, 12])
text(-1.4, 11.7, '(g)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_3*10^-6, el_pitch_3*10^3))

nexttile;
imagesc(x*10^3, z*10^3, dB_TFM_4, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords_4(:,1)*10^3, probe_coords_4(:,2)*10^3, 'wo')
xlim([-1.5, 1.5])
ylim([9, 12])
text(-1.4, 11.7, '(h)', 'Color', 'w')
% title(sprintf('f = %2.1fMHz, p = %3.2fmm', centre_freq_4*10^-6, el_pitch_4*10^3))

xlabel(t2, 'x (mm)')
ylabel(t2, 'z (mm)')

cb = colorbar();
cb.Label.String = 'dB';
cb.Layout.Tile = 'east';