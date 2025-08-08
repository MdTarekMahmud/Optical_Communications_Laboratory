% Graded Index Fiber: Acceptance Angle & Waveguide Dispersion
clc;
clear;
close all;

% --- Fiber parameters (you can change in lab) ---
n1 = 1.48;          % PEAK refractive index at the core center
n2 = 1.46;          % Cladding refractive index
a = 25e-6;          % Core radius (meters)
g = 2;              % Profile parameter (g=2 is a parabolic profile)
lambda = linspace(1.2e-6, 1.6e-6, 200); % Wavelength range (meters)
c = 3e8;            % Speed of light (m/s)

% Part (a): Acceptance Angle
% For a GRIN fiber, this calculates the MAXIMUM acceptance angle,
% which occurs at the center of the core (r=0).
NA = sqrt(n1^2 - n2^2);
theta_a = asind(NA); % in degrees
fprintf('Maximum Numerical Aperture (NA) = %.4f\n', NA);
fprintf('Maximum Acceptance Angle (theta_a) = %.4f degrees\n', theta_a);

% Part (b): Waveguide Dispersion Calculation

% Normalized frequency V
V = (2*pi*a./lambda) .* NA;

% --- Approximation for normalized propagation constant b ---
% For a parabolic profile (g=2) GRIN fiber, a standard approximation is used.
% This is the KEY difference from the step-index code.
b = 1 - (2 ./ V);
b(b < 0) = 0; % Avoid negative values for b

% Effective refractive index
% Using the linear approximation neff ? n2 + b*(n1-n2), valid for small NA
neff = n2 + b.*(n1 - n2);

% Propagation constant beta
beta = (2*pi./lambda) .* neff;

% Numerical 2nd derivative of beta with respect to wavelength
d2beta = gradient(gradient(beta, lambda), lambda);

% Waveguide dispersion (base units: s/m^2)
Dw = -(lambda./c) .* d2beta; % s/m^2

% CORRECT conversion to ps/(nm·km)
Dw_ps = Dw * 1e6; % s/m^2 -> ps/(nm·km)

% Display waveguide dispersion at a specific wavelength (e.g., 1.55 µm)
lambda_given = 1.55e-6;
[~, idx] = min(abs(lambda - lambda_given));
fprintf('Waveguide Dispersion at %.2f µm = %.6f ps/(nm·km)\n', ...
    lambda_given*1e6, Dw_ps(idx));

% Plot
figure;
plot(lambda*1e6, Dw_ps, 'r', 'LineWidth', 1.5); % Changed color to red
xlabel('Wavelength (\mum)');
ylabel('Waveguide Dispersion (ps/(nm·km))');
title('Waveguide Dispersion vs Wavelength (Graded Index Fiber)');
grid on;