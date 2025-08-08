clc;
clear;
close all;
% --- Define constants and Sellmeier coefficients for pure silica ---
lambda_m = linspace(0.3e-6, 2e-6, 1000); % Wavelength range in METERS for calculation
B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
C1_m = (0.0684043e-6)^2; % C coefficients in meters^2
C2_m = (0.1162414e-6)^2;
C3_m = (9.896161e-6)^2;
c = 2.998e8; % Speed of light in m/s
% --- (a) Graph for refractive index (n) vs wavelength (?) ---
% Sellmeier equation using wavelength in meters
n = sqrt(1 + (B1*lambda_m.^2)./(lambda_m.^2 - C1_m) + (B2*lambda_m.^2)./(lambda_m.^2 - C2_m) + (B3*lambda_m.^2)./(lambda_m.^2 - C3_m));

figure;
subplot(2,1,1);
plot(lambda_m * 1e6, n, 'LineWidth', 2); % Plot with lambda in µm
xlim([0.3 2]);
xlabel('Wavelength (\mum)');
ylabel('Refractive Index (n)');
title('a) Refractive Index vs. Wavelength');
grid on;
% --- (b) Graph for material dispersion (Dm) vs wavelength ---
% Calculate the second derivative d²n/d?² numerically
d2n_dlambda2 = gradient(gradient(n, lambda_m), lambda_m);

% CORRECT formula for material dispersion
Dm_si = -(lambda_m./c) .* d2n_dlambda2; % Result is in SI units (s/m^2)

% Convert to standard units of ps/(nm·km)
Dm_ps = Dm_si * 1e6;

subplot(2,1,2);
plot(lambda_m * 1e6, Dm_ps, 'LineWidth', 2, 'Color','r'); % Plot with lambda in µm
xlim([0.3 2]);
ylim([-200 100]); % Set y-axis limits for better visualization
xlabel('Wavelength (\mum)');
ylabel('Dispersion (ps/nm·km)');
title('b) Material Dispersion vs. Wavelength');
grid on;
line([0.3 2], [0 0], 'Color', 'k', 'LineStyle', '--'); % Add zero-dispersion line