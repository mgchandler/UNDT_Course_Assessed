%% UNDT - Exercise 5
% The purpose of this exercise is to design an ultrasonic phased array for
% a specific application and simulate and example image obtained from it.
% This will require the use of various modelling techniques you have
% learned during the course.

clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window

%% Plotting parameters
n = 3;
m = 3;
subfig_label = 'abcdefghijklmnopqrstuvwxyz';
kk = 1;
t = tiledlayout(m,n);

%% Inputs
centre_frequency = [1e6, 2e6, 5e6];
el_pitch = [.1875e-3, .375e-3, .75e-3];
el_separation = .05e-3;
num_els = [64, 32, 16];
velocity_L = 1500;
backwall_distance = 25e-3;
grid_pts = 351;

focal_pt_x = 0.0e-3;
focal_pt_z = 15.00e-3;
          
p = zeros(n*m, grid_pts, grid_pts);

% Generate grid.
x = linspace(-backwall_distance/2, backwall_distance/2, grid_pts);
z = linspace(0, backwall_distance, grid_pts);
dx = x(2) - x(1);
dz = z(2) - z(1);
[X, Z] = meshgrid(x, z);

% Get focal_pt indices
focal_pt_ii = round(focal_pt_x / dx) + round(grid_pts / 2);
focal_pt_kk = round(focal_pt_z / dz);

for ii = 1:n
    for jj = 1:m

        %% Parameters
        wavelength = velocity_L / centre_frequency(ii);
        k = 2 * pi / wavelength;
        omega = 2 * pi * centre_frequency(ii);
        el_width = el_pitch(jj) - el_separation;

        % Mid-coordinates of each element.
        source_x_positions = linspace(0, el_pitch(jj)*(num_els(jj)-1), num_els(jj));
        source_x_positions = source_x_positions - mean(source_x_positions);

        % Adjust X, Z, source_x_positions for vectorisation.
        X_ = repmat(X, 1, 1, length(source_x_positions));
        Z_ = repmat(Z, 1, 1, length(source_x_positions));
        source_positions = reshape(source_x_positions, 1, 1, length(source_x_positions));

        %% Beam Profile of linear array

        % Calculate distance from source to point.
        r_j = sqrt((X_ - source_positions).^2 + (Z_ - .01e-3).^2);
        % Calculate angles made from element to pixel location.
        phi = acos(Z_ ./ r_j);
        directivity_f = el_width * sinc(k * el_width / (2 * pi) * sin(phi));
        % Array time delays. Get the list of distances to this point for all els.
        d_j = squeeze(r_j(focal_pt_kk, focal_pt_ii, :));
        t_j = (d_j - d_j(round(num_els(jj)/2))) / velocity_L;
        % Weight-and-phase delay.
        B_j = repmat(reshape((exp(- 1i * omega * t_j)), 1, 1, num_els(jj)), length(z), length(x), 1);

        % Calculate pressure field. Sum over elements.
        p(kk, :, :) = sum( ...
            1./sqrt(r_j) .* exp(1i*(k * r_j - omega * 0)) .* directivity_f .* B_j, ...
            3 ...
        );
        kk = kk + 1;
    end
end

kk = 1;

for ii = 1:n
    for jj = 1:m
        % Get aspect ratios.
        this_p = abs(reshape(p(kk, 11:end, :), grid_pts-10, grid_pts));
        max_ = max(this_p, [], 'all');
        dB_p = 20 * log10(this_p ./ max_);
        is_focus = logical(dB_p > -3);
        
        is_x = x(any(is_focus, 1))*10^3;
        is_z = z(any(is_focus, 2))*10^3;
        
        aspect_ratio = (is_x(end) - is_x(1)) / (is_z(end) - is_z(1));
        
        source_x_positions = linspace(0, el_pitch(jj)*(num_els(jj)-1), num_els(jj));
        source_x_positions = (source_x_positions - mean(source_x_positions)) * 10^3;
        nexttile;
        hold on
        box on
        imagesc(x*10^3, z*10^3, abs(reshape(p(kk, :, :), grid_pts, grid_pts)))
        scatter(source_x_positions, zeros(num_els(jj),1), 'wo')
        scatter(focal_pt_x*10^3, focal_pt_z*10^3, 'r.')
        title(sprintf('f = %2.1fMHz\np = %3.2fmm, n = %d', centre_frequency(ii)*10^-6, el_pitch(jj)*10^3, num_els(jj)))
        text(-12.5, 21, sprintf('(%s) λ = %2.1fmm\nAR = %3.2f', subfig_label(kk+9), 10^3*velocity_L / centre_frequency(ii), aspect_ratio), 'Color', 'white')
        xlim([x(1)*10^3, x(end)*10^3])
        ylim([z(1)*10^3, z(end)*10^3])
        kk = kk + 1;
    end
end
xlabel(t, 'x (mm)')
ylabel(t, 'z (mm)')



%% Plotting parameters
n = 5;
m = 2;
subfig_label = 'abcdefghijklmnopqrstuvwxyz';
kk = 1;
figure(2)
t = tiledlayout(m,n);

%% Inputs
centre_frequency = 5e6;
el_pitch = .15e-3;
el_separation = .05e-3;
num_els = 64;
velocity_L = 1500;
backwall_distance = 20e-3;
grid_pts = 351;

focal_pt_x = [-2.25e-3, ...
              -1.00e-3, ...
               0.25e-3, ...
               1.50e-3, ...
               2.75e-3, ...
               0.0, ...
               0.0, ...
               0.0, ...
               0.0, ...
               0.0];
focal_pt_z = [ 5.00e-3, ...
               7.50e-3, ...
              10.00e-3, ...
              12.50e-3, ...
              15.00e-3, ...
              5.00e-3, ...
               7.50e-3, ...
              10.00e-3, ...
              12.50e-3, ...
              15.00e-3];
          
p = zeros(n*m, grid_pts, grid_pts);

for ii = 1:n*m

    %% Parameters
    wavelength = velocity_L / centre_frequency;
    k = 2 * pi / wavelength;
    omega = 2 * pi * centre_frequency;
    el_width = el_pitch - el_separation;
    % Generate grid.
    x = linspace(-backwall_distance/2, backwall_distance/2, grid_pts);
    z = linspace(0, backwall_distance, grid_pts);
    dx = x(2) - x(1);
    dz = z(2) - z(1);
    [X, Z] = meshgrid(x, z);

    % Get focal_pt indices
    focal_pt_ii = round(focal_pt_x(ii) / dx) + round(grid_pts / 2);
    focal_pt_kk = round(focal_pt_z(ii) / dz);

    % Mid-coordinates of each element.
    source_x_positions = linspace(0, el_pitch*(num_els-1), num_els);
    source_x_positions = source_x_positions - mean(source_x_positions);

    % Adjust X, Z, source_x_positions for vectorisation.
    X = repmat(X, 1, 1, length(source_x_positions));
    Z = repmat(Z, 1, 1, length(source_x_positions));
    source_positions = reshape(source_x_positions, 1, 1, length(source_x_positions));

    %% Beam Profile of linear array

    % Calculate distance from source to point.
    r_j = sqrt((X - source_positions).^2 + (Z - .01e-3).^2);
    % Calculate angles made from element to pixel location.
    phi = acos(Z ./ r_j);
    directivity_f = el_width * sinc(k * el_width / (2 * pi) * sin(phi));
    % Array time delays. Get the list of distances to this point for all els.
    d_j = squeeze(r_j(focal_pt_kk, focal_pt_ii, :));
    t_j = (d_j - d_j(round(num_els/2))) / velocity_L;
    % Weight-and-phase delay.
    B_j = repmat(reshape((exp(- 1i * omega * t_j)), 1, 1, num_els), length(z), length(x), 1);

    % Calculate pressure field. Sum over elements.
    p(ii, :, :) = sum( ...
        1./sqrt(r_j) .* exp(1i*(k * r_j - omega * 0)) .* directivity_f .* B_j, ...
        3 ...
    );

end

for ii = 1:length(focal_pt_x)
    nexttile;
    hold on
    box on
    imagesc(x*10^3, z*10^3, abs(reshape(p(ii, :, :), grid_pts, grid_pts)), [min(abs(p), [], 'all'), max(abs(p), [], 'all')])
    scatter(source_x_positions*10^3, zeros(num_els,1), 'wo')
    scatter(focal_pt_x([1:length(focal_pt_x)]~=ii)*10^3, focal_pt_z([1:length(focal_pt_x)]~=ii)*10^3, 'ro')
    scatter(focal_pt_x([1:length(focal_pt_x)]==ii)*10^3, focal_pt_z([1:length(focal_pt_x)]==ii)*10^3, 'bo')
    title(sprintf('x = %3.2fmm\nz = %3.2fmm', focal_pt_x(ii)*10^3, focal_pt_z(ii)*10^3))
    text(10^3*focal_pt_x(ii) - 2.8, 10^3*focal_pt_z(ii) + 2.7, sprintf('(%s)', subfig_label(ii)), 'Color', 'white')
    xlim([10^3*focal_pt_x(ii) - 3, 10^3*focal_pt_x(ii) + 3])
    ylim([10^3*focal_pt_z(ii) - 3, 10^3*focal_pt_z(ii) + 3])
end

t.XLabel.String = 'x (mm)';
t.YLabel.String = 'z (mm)';
t.TileSpacing = 'tight';
cb = colorbar;
cb.Layout.Tile = 'east';