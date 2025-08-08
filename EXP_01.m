clc;
clear;
close all;
% --- Fiber Parameters ---
n1 = 1.48;          % Core refractive index
n2 = 1.46;          % Cladding refractive index
a = 4.1e-6;         % Core radius (in meters)
lambda = linspace(0.8e-6, 1.6e-6, 500);  % Wavelength range in meters
c = 3e8;            % Speed of light (m/s)
% --- (a) Calculations ---
% Numerical Aperture (NA)
NA = sqrt(n1^2 - n2^2);
% V-number for each wavelength
V = (2 * pi * a ./ lambda) * NA;
% --- CORRECTED Propagation constant ? calculation ---
% 1. First, find the normalized propagation constant 'b' for the LP01 mode.
%    The Gloge formula is a standard and accurate approximation.
b = (1.1428 - 0.996./V).^2;
b(V <= 1.5) = 0; % Gloge formula is valid for V > 1.5
% 2. Now, calculate the correct, wavelength-dependent 'neff' using b.
neff = sqrt(n2^2 + b .* (n1^2 - n2^2));
% 3. Finally, calculate the correct 'beta'.
beta = (2 * pi ./ lambda) .* neff;
% --- End of correction ---
% Cutoff V-number for the fundamental mode (LP11 cutoff)
V_cutoff = 2.405;
% Cutoff wavelength (?c) for single mode operation
lambda_c = (2 * pi * a * NA) / V_cutoff;
% Display Results
fprintf('Numerical Aperture (NA): %.4f\n', NA);
fprintf('Cutoff Wavelength (?c): %.4e meters (%.1f nm)\n\n', lambda_c, lambda_c*1e9);
% Display Table of Sample Values
fprintf('  Wavelength(nm)      V-number      Propagation Constant ? (rad/m)\n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:50:length(lambda) % Print every 50th value to limit output
    fprintf('%10.2f           %8.4f         %15.4e\n', lambda(i)*1e9, V(i), beta(i));
end
% --- (b) Plot: Propagation Constant (?) vs V-number ---
figure;
plot(V, beta, 'LineWidth', 2);
xlabel('V-number');
ylabel('Propagation constant ? (rad/m)');
title('Propagation Constant ? vs V-number');
grid on;
% --- (c) Plot: Wavelength (?) vs V-number ---
figure;
plot(lambda * 1e9, V, 'LineWidth', 2, 'Color', 'r');  % ? in nm
xlabel('Wavelength (nm)');
ylabel('V-number');
title('Wavelength vs V-number');
grid on;