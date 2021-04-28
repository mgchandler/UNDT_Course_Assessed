%% UNDT Exercise 5b
% Part 3

%% Inputs

centre_frequency = 2.5e6; % Probe frequency
scatterer_coords = [[-2.25e-3,  5.00e-3]; ...
                    [-1.00e-3,  7.50e-3]; ...
                    [ 0.25e-3, 10.00e-3]; ...
                    [ 1.50e-3, 12.50e-3]; ...
                    [ 2.75e-3, 15.00e-3]; ...
                    [ 0.00e-3,  5.00e-3]; ...
                    [ 0.00e-3,  7.50e-3]; ...
                    [ 0.00e-3, 10.00e-3]; ...
                    [ 0.00e-3, 12.50e-3]; ...
                    [ 0.00e-3, 15.00e-3]];
num_els = 32;
el_pitch = .1e-3;
backwall_depth = 20.0e-3;
grid_pts = 201;
velocity_L = 1500.0;
num_cycles = 5;
          


%% Parameters
probe_coords = linspace(0, el_pitch*(num_els-1), num_els);
probe_coords = probe_coords - mean(probe_coords);

omega = 2 * pi * centre_frequency;

num_scatterers = size(scatterer_coords, 1);

x = linspace(-backwall_depth/2, backwall_depth/2, grid_pts);
z = linspace(0, backwall_depth, grid_pts);
[X, Z] = meshgrid(x, z);

max_t = 2.1 * sqrt((x(end) - x(1))^2 + (z(end) - z(1))^2) / velocity_L;
dt = 1 / (centre_frequency * 5);
time = [0 : dt : max_t];

half_pulse = 5 / (centre_frequency * 2);

[tx, rx] = meshgrid([1:num_els]);
tx = reshape(tx, length(tx)^2, 1);
rx = reshape(rx, length(rx)^2, 1);



%% Input spectrum
in_sig = sin(2 * pi * centre_frequency * time)' .* ...
         fn_hanning(length(time), half_pulse / max_t, half_pulse / max_t);
fft_pts = 2^nextpow2(length(time));

in_spec = fft(in_sig, fft_pts);
in_spec = in_spec(1:round(fft_pts/2));

df = 1/(time(end) - time(1));
freq = [0 : df : df*(length(in_spec)-1)];



%% Ray tracing
sc_dists = (sqrt((probe_coords(tx)' - scatterer_coords(:, 1)').^2 + ...
                 (0.01e-3 - scatterer_coords(:, 2)').^2) + ...
            sqrt((probe_coords(rx)' - scatterer_coords(:, 1)').^2 + ...
                 (0.01e-3 - scatterer_coords(:, 2)').^2));
sc_times = sc_dists / velocity_L;
bw_dists = (sqrt((probe_coords(tx)' - probe_coords(rx)').^2 + ...
                 (2*(backwall_depth + 0.01e-3)).^2));
bw_times = bw_dists / velocity_L;
       
       
       
%% Propagate signal
out_spec = 0;
for scatterer = 1:num_scatterers
    diag1 = spdiags(in_spec, 0, length(in_spec), length(in_spec));
    delta = exp(-1i * 2*pi*freq' * sc_times(:, scatterer)') * 0.01;
    diag2 = spdiags(sqrt(sc_dists(:, scatterer)), 0, length(sc_dists(:, scatterer)), length(sc_dists(:, scatterer)));
    out_spec = out_spec + ( ...
        spdiags(in_spec, 0, length(in_spec), length(in_spec)) * ...
        exp(-1i * 2*pi*freq' * sc_times(:, scatterer)') * 0.01 / ...
        spdiags(sqrt(sc_dists(:, scatterer)), 0, length(sc_dists(:, scatterer)), length(sc_dists(:, scatterer))) ...
    );
end
% out_spec = out_spec + ( ...
%     in_spec' .* exp(-1i * omega * bw_times) ./ sqrt(bw_dists)...
% );

%% Convert to timetraces
FMC_timetraces = ifft(out_spec, fft_pts, 1);
FMC_timetraces = FMC_timetraces(1:length(time), :);

figure(1)
imagesc([1:num_els^2], time, abs(FMC_timetraces))

%% TFM

lookup_times = (sqrt((probe_coords(tx)' - X(:)').^2 + ...
                     (0.01e-3 - Z(:)').^2) + ...
                sqrt((probe_coords(rx)' - X(:)').^2 + ...
                     (0.01e-3 - Z(:)').^2) ...
    ) / velocity_L;
Im = 0;
for tr = 1:num_els^2
    tau = reshape(lookup_times(tr, :), size(X, 1), size(X, 2));
    Im = Im + ...
        interp1(time, FMC_timetraces(:, tr), tau, 'linear', 0);
end

max_ = max(abs(Im), [], 'all');
Im = 20 * log10(abs(Im) ./ max_);

figure(2)
imagesc(x, z, Im)
hold on
scatter(scatterer_coords(:, 1), scatterer_coords(:, 2), 'ro')