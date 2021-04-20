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

signals = ifft(hanning_spectra, fft_pts, 1);
signals = signals(1:length(time), :);

% Work out response locations.

threshold = min(max(abs(signals), [], 1)) / 175;
is_response = logical(abs(signals) > threshold);

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
    if and(start_idxs(ii, 1) ~= 0, start_idxs(ii, 2) ~= 0)
        bearing_thickness(ii) = (time(start_idxs(ii, 2)) - time(start_idxs(ii, 1))) * 5900.0 / 2;
    end
end

figure(1)
c = get(gca, 'colororder');
hold on
shift = 0;
for ii = 1:length(pos)
    plot(time, abs(signals(:, ii)) + shift, 'Color', c(mod(ii-1, 7)+1, :));
    plot([time(1), time(end)], [threshold + shift, threshold + shift], 'Color', c(mod(ii-1, 7)+1, :))
    shift = shift + .1;
end

figure(2)
scatter(pos, bearing_thickness)

