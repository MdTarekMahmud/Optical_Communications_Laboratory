clc;
clear;
close all;
% --- Fiber Parameters for Graded Index Fiber ---
n1 = 1.48;          % PEAK core refractive index at the center
n2 = 1.46;          % Cladding refractive index
a = 25e-6;          % Core radius (in meters)
g = 2;              % Profile parameter (2 for parabolic profile)
lambda = linspace(0.8e-6, 1.6e-6, 500);  % Wavelength range (in meters)
% --- (a) Calculations ---
% Numerical Aperture (NA)
NA = sqrt(n1^2 - n2^2);
% V-number
V = (2 * pi * a ./ lambda) .* NA;
% --- CORRECTED Propagation constant ? calculation ---
% 1. Use a standard approximation for the normalized propagation constant 'b'
%    for a parabolic (g=2) GRIN fiber.
b = 1 - (2 ./ V);
b(b < 0) = 0; % b cannot be negative

% 2. Calculate the correct, wavelength-dependent 'neff' from b.
%    The approximation neff ? n2 + b*(n1-n2) is accurate for small NA.
neff = n2 + b .* (n1 - n2);

% 3. Calculate the correct 'beta'
beta = (2 * pi ./ lambda) .* neff;
% --- End of correction ---

% --- CORRECTED Cutoff V-number for a parabolic (g=2) GRIN fiber ---
V_cutoff = 3.518; 

% Cutoff wavelength ?c
lambda_c = (2 * pi * a * NA) / V_cutoff;

% Display Results
fprintf('Graded Index Fiber Analysis:\n');
fprintf('Numerical Aperture (NA): %.4f\n', NA);
fprintf('Cutoff Wavelength (?c): %.4e meters (%.1f nm)\n\n', lambda_c, lambda_c*1e9);

% Display Sample Data
fprintf('  Wavelength(nm)      V-number      Propagation Constant ? (rad/m)\n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:50:length(lambda)  % Print every 50th value
    fprintf('%10.2f           %8.4f         %15.4e\n', lambda(i)*1e9, V(i), beta(i));
end
% --- (b) Plot: Propagation Constant vs V-number ---
figure;
plot(V, beta, 'LineWidth', 2);
xlabel('V-number');
ylabel('Propagation Constant ? (rad/m)');
title('Propagation Constant vs V-number (Graded Index Fiber)');
grid on;
% --- (c) Plot: Wavelength vs V-number ---
figure;
plot(lambda*1e9, V, 'LineWidth', 2, 'Color', 'r');
xlabel('Wavelength (nm)');
ylabel('V-number');
title('Wavelength vs V-number (Graded Index Fiber)');
grid on;