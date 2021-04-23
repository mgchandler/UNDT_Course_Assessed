clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window

load('bearing_casing_bscan.mat')

% Filter data and get envelope.

fft_pts = 2^nextpow2(length(time));
spectra = fft(voltage, fft_pts, 1);
spectra = spectra(1:fft_pts/2, :);

df = 1/(time(end) - time(1));
freq = [0 : df : df*(fft_pts/2 - 1)];

gaussian_window = fn_gaussian(fft_pts/2, 5e6/freq(end), 5e6/freq(end));
gaussian_spectra = spectra .* gaussian_window;
hanning_window = fn_hanning(fft_pts/2, 5e6/freq(end), 5e6/freq(end));
hanning_spectra = spectra .* hanning_window;
hanning_hi_window = fn_hanning_hi_pass(length(spectra), 2*5e6/freq(end), 3*5e6/freq(end));
hanning_hi_spectra = spectra .* hanning_hi_window;
hanning_lo_window = fn_hanning_lo_pass(length(spectra), 2*5e6/freq(end), 3*5e6/freq(end));
hanning_lo_spectra = spectra .* hanning_lo_window;

g_signals = ifft(gaussian_spectra, fft_pts, 1);
g_signals = g_signals(1:length(time), :);

h_signals = ifft(hanning_spectra, fft_pts, 1);
h_signals = h_signals(1:length(time), :);

hh_signals = ifft(hanning_hi_spectra, fft_pts, 1);
hh_signals = hh_signals(1:length(time), :);

hl_signals = ifft(hanning_lo_spectra, fft_pts, 1);
hl_signals = hl_signals(1:length(time), :);

% Plot filtered data.

figure(1)
subplot(2,2,1)
imagesc(time*10^6, pos*10^3, abs(g_signals'))
xlabel('Time (\mus)')
ylabel('Position (mm)')

subplot(2,2,2)
imagesc(time*10^6, pos*10^3, abs(h_signals'))
xlabel('Time (\mus)')
ylabel('Position (mm)')

subplot(2,2,3)
imagesc(time*10^6, pos*10^3, abs(hh_signals'))
xlabel('Time (\mus)')
ylabel('Position (mm)')

subplot(2,2,4)
imagesc(time*10^6, pos*10^3, abs(hl_signals'))
xlabel('Time (\mus)')
ylabel('Position (mm)')

c = colorbar();
c.Label.String = 'Voltage |V(x, t)|';

% Work out response locations.

threshold = min(max(abs(h_signals), [], 1)) / 200;
is_response = logical(abs(h_signals) > threshold);

lengths = zeros(length(pos), 1);
start_idxs = zeros(length(pos), 3);

for ii = 1:length(pos)
    end_idxs = [0];
    this_response_start = 1;
    this_response_end = 2;
    for jj = 1:length(time)-1
        if and(and(~is_response(jj, ii), is_response(jj+1, ii)), jj > end_idxs(end)+100)
            % If at the leading edge
            start_idxs(ii, this_response_start) = jj;
            this_response_start = this_response_start + 1;
        end
        if and(is_response(jj, ii), ~is_response(jj+1, ii))
            % If at the trailing edge
            end_idxs(this_response_end) = jj+1;
            this_response_end = this_response_end + 1;
        end
    end
end

% Calculate thickness from response locations.

bearing_thickness = zeros(length(pos), 1);
for ii = 1:length(pos)
    % Confirm that there are at least two leading edges of responses.
    if and(start_idxs(ii, 1) ~= 0, start_idxs(ii, 2) ~= 0)
        bearing_thickness(ii) = (time(start_idxs(ii, 2)) - time(start_idxs(ii, 1))) * 5900.0 / 2;
    end
end
% Remove any data where there are not two leading edges.
pos_ = pos(bearing_thickness ~= 0);
bearing_thickness(bearing_thickness == 0) = [];
h_signals_ = h_signals(:, bearing_thickness ~= 0);

thickness_at_125 = interp1(pos_, bearing_thickness, 0.125, 'linear');

% Plot thickness against transducer position.

figure(2)
scatter(pos_*10^3, bearing_thickness*10^3, 'r.')
hold on
scatter([125], [thickness_at_125*10^3], 'bx')
text(125*1.025, thickness_at_125*10^3*0.99, sprintf('%4.1fmm', thickness_at_125*10^3))
xlabel('Transducer Position (mm)')
ylabel('Bearing Thickness (mm)')
box on