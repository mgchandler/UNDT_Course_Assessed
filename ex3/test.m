clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window
 
%% Inputs
 
num_els = 64;
centre_freq = 5e6;
el_pitch = .15e-3;
el_sep = .05e-3;
v_L = 1500.0;
backwall_dist = 20.0e-3;
grid_pts = 251;
imaging_aperture = el_pitch * num_els / 4;
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
 
lambda = v_L / centre_freq;
 
probe_coords = zeros(num_els, 2);
probe_coords(:, 1) = linspace(0, el_pitch*(num_els-1), num_els);
probe_coords(:, 1) = probe_coords(:, 1) - mean(probe_coords(:, 1));
 
[tx, rx] = meshgrid([1:num_els]);
tx = reshape(tx, length(tx)^2, 1);
rx = reshape(rx, length(rx)^2, 1);
 
%% Get input time signal and freq spectrum
 
max_time = 2.1 * sqrt(2 * backwall_dist^2) / v_L;
dt = 1 / centre_freq / oversample;
time = [0 : dt : max_time];
 
half_pulse = 5 / (2 * centre_freq);
 
input_signal = sin(2 * pi * centre_freq * time)' .* ...
               fn_hanning(length(time), half_pulse / max_time, half_pulse / max_time);
           
% Bring peak of input pulse to time = 0.
time = time - half_pulse;
           
           
           
fft_pts = 2^nextpow2(length(time));
 
df = 1 / max_time;
freq = [0 : df : df*(fft_pts/2-1)];
omega = 2 * pi * freq;
 
input_spectrum = fft(input_signal, fft_pts);
input_spectrum = input_spectrum(1:fft_pts/2);
 
%% Do ray tracing
 
scat_dists_tx = sqrt((probe_coords(tx, 1) - scat_coords(:, 1)').^2 + ...
                      (probe_coords(tx, 2) - scat_coords(:, 2)').^2);
scat_dists_rx = sqrt((probe_coords(rx, 1) - scat_coords(:, 1)').^2 + ...
                     (probe_coords(rx, 2) - scat_coords(:, 2)').^2);
scat_times = (scat_dists_tx + scat_dists_rx) / v_L;
                      
bw_dists = (sqrt((probe_coords(tx, 1)' - probe_coords(rx, 1)').^2 + ...
                 (2*(backwall_dist - probe_coords(tx, 2)')).^2));
bw_times = bw_dists / v_L;
 
% Get angles for directivity.
scat_theta_tx = acos((scat_coords(:, 2)' - probe_coords(tx, 2)) ./ ...
                     sqrt((probe_coords(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords(tx, 2) - scat_coords(:, 2)').^2));
scat_theta_rx = acos((scat_coords(:, 2)' - probe_coords(rx, 2)) ./ ...
                     sqrt((probe_coords(tx, 1) - scat_coords(:, 1)').^2 + ...
                          (probe_coords(tx, 2) - scat_coords(:, 2)').^2));
                      
bw_theta_tx = acos((2*(backwall_dist - probe_coords(tx, 2)')) ./ bw_dists);
bw_theta_rx = acos((2*(backwall_dist - probe_coords(rx, 2)')) ./ bw_dists);
 
%% Propagate signal
 
% Get scatterer and backwall signal amplitudes.
 
% Directivity includes both transmit and receive paths. N.B. this function
% takes into account that in Matlab, sinc(x) := sin(πx)/( πx)
scat_dir = (el_pitch - el_sep)^2 * ...
           sinc((el_pitch - el_sep) / lambda * sin(scat_theta_tx)) .* ...
           sinc((el_pitch - el_sep) / lambda * sin(scat_theta_rx));
% Amplitude includes scattering, directivity and beam spreading.
scat_amp = 0.01 * scat_dir ./ (sqrt(scat_dists_tx) .* sqrt(scat_dists_rx));
                               
bw_dir = (el_pitch - el_sep)^2 * ...
         sinc((el_pitch - el_sep) / lambda * sin(bw_theta_tx)) .* ...
         sinc((el_pitch - el_sep) / lambda * sin(bw_theta_rx));
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
 
FMC_data = ifft(output_spectra, fft_pts, 1);
FMC_data = FMC_data(1:length(time), :); 
FMC_dt = 1 / (fft_pts * abs(freq(2) - freq(1)));
FMC_time = [0 : FMC_dt : FMC_dt * (length(time)-1)];
 
% End of simulation of wave.
 
%% Imaging - TFM
 
x = linspace(- backwall_dist / 2 - 1e-3, backwall_dist / 2 + 1e-3, grid_pts);
z = linspace(-1e-3, backwall_dist + 1e-3, grid_pts);
[X, Z] = meshgrid(x, z);
 
%% Plane
 
Plane_lookup_times = 2 * Z(:)' / v_L;
% Work out if each point is below the aperture for each tr pair
is_in_aperture = logical( ...
    abs(probe_coords(:, 1) - X(:)') <= imaging_aperture / 2 ...
);
 
Im_Plane = 0;
tr = 1;
for ti = 1:num_els
    for ri = 1:num_els
        % Skip all tr pairs where there are no contributions.
        if and(any(is_in_aperture(ti, :)), any(is_in_aperture(ri, :)))
            Im_Plane = Im_Plane + reshape( ...
    ... % Only contribute if this point in the image is below both tx and rx.
                is_in_aperture(ti, :) .* is_in_aperture(ri, :) .* ...
    ... % Interpolate to find the data.
                interp1(FMC_time, FMC_data(:, tr), Plane_lookup_times, 'linear', 0), ...
                size(X, 1), size(X, 2) ...
            );
        end
        tr = tr + 1;
    end
end
 
max_ = max(abs(Im_Plane), [], 'all');
dB_Plane = 20 * log10(abs(Im_Plane) ./ max_);
 
%% Focussed
 
Focus_lookup_times = (sqrt((probe_coords(tx, 1) - X(:)').^2 + ...
                           (probe_coords(tx, 2) - Z(:)').^2) + ...
                      sqrt((probe_coords(rx, 1) - X(:)').^2 + ...
                           (probe_coords(rx, 2) - Z(:)').^2) ...
    ) / v_L;
 
Im_Focus = 0;
tr = 1;
for ti = 1:num_els
    for ri = 1:num_els
        if and(any(is_in_aperture(ti, :)), any(is_in_aperture(ri, :)))
            Im_Focus = Im_Focus + reshape( ...
                is_in_aperture(ti, :) .* is_in_aperture(ri, :) .* ...
                interp1(FMC_time, FMC_data(:, tr), Focus_lookup_times(tr, :), 'linear', 0), ...
                size(X, 1), size(X, 2));
        end
        tr = tr + 1;
    end
end
 
max_ = max(abs(Im_Focus), [], 'all');
dB_Focus = 20 * log10(abs(Im_Focus) ./ max_);
 
%% Sector
 
r = linspace(0, 1.1 * sqrt(5)/2 * backwall_dist, grid_pts);
theta = linspace(-pi/4, pi/4, grid_pts);
[R, Theta] = meshgrid(r, theta);
 
Sector_lookup_times = ( ...
    2*R(:)' + probe_coords(tx, 1) * sin(Theta(:)') + probe_coords(rx, 1) * sin(Theta(:)') ...
) / v_L;
tau_Sector = reshape(Sector_lookup_times, num_els^2, size(R, 1), size(R, 2));
 
Im_Sector = 0;
for tr = 1:num_els^2
    Im_Sector = Im_Sector + ...
        interp1(FMC_time, FMC_data(:, tr), squeeze(tau_Sector(tr, :, :)), 'linear', 0);
end
 
max_ = max(abs(Im_Sector), [], 'all');
dB_Sector = 20 * log10(abs(Im_Sector) ./ max_);
dB_Sector(dB_Sector < -40) = -40;
 
%% TFM
 
% TFM lookup times are identical to focussed lookup times. Difference is
% that we focus everywhere, instead of just below the aperture.
 
Im_TFM = 0;
tau_TFM = reshape(Focus_lookup_times, num_els^2, size(X, 1), size(X, 2));
for tr = 1:num_els^2
    Im_TFM = Im_TFM + ...
        interp1(FMC_time, FMC_data(:, tr), squeeze(tau_TFM(tr, :, :)), 'linear', 0);
end
 
max_ = max(abs(Im_TFM), [], 'all');
dB_TFM = 20 * log10(abs(Im_TFM) ./ max_);
 
%% Plotting
 
figure(1)
subplot(1,2,1)
xlim([min(x)*10^3, max(x)*10^3])
ylim([min(z)*10^3, max(z)*10^3])
hold on
plot([-backwall_dist*10^3/2, backwall_dist*10^3/2], [backwall_dist*10^3, backwall_dist*10^3], 'b')
scatter(probe_coords(:, 1)*10^3, probe_coords(:, 2)*10^3, 'ko')
scatter(scat_coords(:, 1)*10^3, scat_coords(:, 2)*10^3, 'r.')
xlabel('x (mm)')
ylabel('y (mm)')
box on
text(-10, 19, '(a)')
subplot(1,2,2)
imagesc(time*10^6, [1:num_els^2], abs(FMC_data'))
xlabel('Time (μs)')
ylabel('Tx-Rx pair')
box on
text(1, 120, '(b)', 'Color', 'w')
cb = colorbar();
cb.Label.String = 'Signal';
 
figure(2)
tiledlayout(2, 2)
nexttile;
imagesc(x*10^3, z*10^3, dB_Plane, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords(:,1)*10^3, probe_coords(:,2)*10^3, 'wo')
xlabel('x (mm)')
ylabel('z (mm)')
text(-10, 19, '(a)', 'Color', 'w')
title('Plane Scan')
 
nexttile;
imagesc(x*10^3, z*10^3, dB_Focus, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords(:,1)*10^3, probe_coords(:,2)*10^3, 'wo')
xlabel('x (mm)')
ylabel('z (mm)')
text(-10, 19, '(b)', 'Color', 'w')
title('Focussed Scan')
 
nexttile;
dBScaled = round((dB_Sector(:)+40));
map = colormap(parula(41));
polarscatter(Theta(:), R(:)*10^3, [], map(dBScaled+1,:), 'filled')
thetalim([-90, 90])
ax = gca;
ax.ThetaAxisUnits = 'radians';
ax.ThetaZeroLocation = 'top';
ax.RAxis.Label.String = 'r (mm)';
text(pi/4, 23, '(c)', 'Color', 'w')
title('Sector Scan')
 
nexttile;
imagesc(x*10^3, z*10^3, dB_TFM, [-40, 0])
ax = gca();
set(ax, 'YDir', 'normal')
hold on
scatter(probe_coords(:,1)*10^3, probe_coords(:,2)*10^3, 'wo')
xlabel('x (mm)')
ylabel('z (mm)')
text(-10, 19, '(d)', 'Color', 'w')
title('TFM')
        % Get aspect ratios. Note: this method only works when focal_pt_x 
        % = 0, and will likely change for different focal_pt_z. Keep z 
        % constant to make comparison valid.
cb = colorbar();
cb.Label.String = 'dB';
cb.Layout.Tile = 'east';

