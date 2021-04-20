%% UNDT - Exercise 5
% The purpose of this exercise is to design an ultrasonic phased array for
% a specific application and simulate and example image obtained from it.
% This will require the use of various modelling techniques you have
% learned during the course.

%% Inputs
centre_frequency = 1e6;
el_pitch = .4e-3;
el_separation = .05e-3;
num_els = 32;
velocity_L = 1500;
backwall_distance = 20e-3;
grid_pts = 201;

focal_pt_x = 0e-3;
focal_pt_z = 10e-3;

%% Parameters
wavelength = velocity_L / centre_frequency;
k = 2 * pi / wavelength;
omega = 2 * pi * centre_frequency;
el_width = el_pitch - el_separation;
% Generate grid.
x = linspace(-backwall_distance / 2, backwall_distance / 2, grid_pts);
z = linspace(0, backwall_distance, grid_pts);
dx = x(2) - x(1);
dz = z(2) - z(1);
[X, Z] = meshgrid(x, z);

% Get focal_pt indices
focal_pt_ii = round(focal_pt_x / dx) + round(grid_pts / 2);
focal_pt_kk = round(focal_pt_z / dz);

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
p = sum( ...
    1./sqrt(r_j) .* exp(1i*(k * r_j - omega * 0)) .* directivity_f .* B_j, ...
    3 ...
);

imagesc(x, z, abs(p))