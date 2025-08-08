clc; clear; close all;

% Wavelength range in micrometers (0.4 µm to 2.0 µm)
lambda = linspace(0.4, 2.0, 500);  % µm

% Constants for loss modeling (approximate typical values)
A = 0.8;              % absorption coefficient constant (empirical)
lambda0 = 1.1;        % reference wavelength for absorption (µm)
B = 0.003;            % Rayleigh scattering constant

% Calculate Absorption Loss (alpha_a)
alpha_a = A * exp(-lambda / lambda0);

% Calculate Rayleigh Scattering Loss (alpha_s)
alpha_s = B ./ (lambda.^4);

% Total Loss
alpha_total = alpha_a + alpha_s;

% Plotting
figure;
plot(lambda, alpha_a, 'b--', 'LineWidth', 1.5); hold on;
plot(lambda, alpha_s, 'r--', 'LineWidth', 1.5);
plot(lambda, alpha_total, 'k-', 'LineWidth', 2);

xlabel('Wavelength (µm)');
ylabel('Loss (dB/km)');
title('Loss vs. Wavelength for Silica Optical Fiber');
legend('Absorption Loss', 'Scattering Loss', 'Total Loss');
grid on;

% Display total loss in the Command Window
fprintf('Wavelength (µm)   Absorption (dB/km)   Scattering (dB/km)   Total Loss (dB/km)\n');
fprintf('----------------  ------------------  -------------------  -------------------\n');

% Print selected values (every 25th point for readability)
for i = 1:25:length(lambda)
    fprintf('%0.3f               %0.4f               %0.4f              %0.4f\n', ...
        lambda(i), alpha_a(i), alpha_s(i), alpha_total(i));
end
