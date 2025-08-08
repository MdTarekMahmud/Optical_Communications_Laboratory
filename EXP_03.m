% Step Index Fiber: Acceptance Angle & Waveguide Dispersion (IMPROVED)
clc;
clear;
close all;

% --- Fiber parameters (you can change in lab) ---
n1 = 1.48;          % Core refractive index
n2 = 1.46;          % Cladding refractive index
a = 4.5e-6;         % Core radius (meters)
lambda = linspace(1.2e-6, 1.6e-6, 200); % Wavelength range (meters)
c = 3e8;            % Speed of light (m/s)

% Part (a): Acceptance Angle (Correct)
NA = sqrt(n1^2 - n2^2);
theta_a = asind(NA); % in degrees
fprintf('Numerical Aperture (NA) = %.4f\n', NA);
fprintf('Acceptance Angle (theta_a) = %.4f degrees\n', theta_a);

% Part (b): Waveguide Dispersion Calculation (Improved Method)

% Normalized frequency V
V = (2*pi*a./lambda) .* NA;

% --- IMPROVED approximation for normalized propagation constant b ---
% Using the Gloge empirical formula, which is more accurate for the LP01 mode
b = (1.1428 - 0.996./V).^2;
b(V <= 1.5) = 0; % Gloge formula is valid for V > 1.5

% --- Use the standard definition of Waveguide Dispersion ---
% Dw = - (lambda/c) * d^2(neff)/d(lambda)^2

% Effective refractive index
neff = sqrt(b.*(n1^2 - n2^2) + n2^2);

% Numerical 2nd derivative of neff with respect to wavelength
d2neff = gradient(gradient(neff, lambda), lambda);

% Waveguide dispersion using the standard formula (base units: s/m^2)
Dw = -(lambda./c) .* d2neff;

% CORRECT conversion to ps/(nm·km)
Dw_ps = Dw * 1e6;  % s/m^2 -> ps/(nm·km)

% Display waveguide dispersion at a specific wavelength (e.g., 1.55 µm)
lambda_given = 1.55e-6;
[~, idx] = min(abs(lambda - lambda_given));
fprintf('Waveguide Dispersion at %.2f µm = %.6f ps/(nm·km)\n', ...
    lambda_given*1e6, Dw_ps(idx));

% Plot
figure;
plot(lambda*1e6, Dw_ps, 'b', 'LineWidth', 1.5);
xlabel('Wavelength (\mum)');
ylabel('Waveguide Dispersion (ps/(nm·km))');
title('Waveguide Dispersion vs Wavelength (Step Index Fiber)');
grid on;