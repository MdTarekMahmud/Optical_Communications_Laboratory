clc;
clear;
close all;

% Fiber Parameters
n1_0 = 1.5;           % Refractive index at core center
n2 = 1.48;            % Cladding refractive index
a = 25e-6;            % Core radius (25 µm)
Delta = (n1_0 - n2) / n1_0;  % Relative index difference

r = linspace(0, a, 500);   % Radial positions from 0 to a

% Profile parameters to test
g_values = [1, 2, 10];      % Triangular, parabolic, near step-index

% Create figure with 2 subplots
figure('Name', 'Graded Index Fiber Profiles');

% Subplot 1: Acceptance Angle
subplot(2,1,1);
hold on;
grid on;
title('a) Acceptance Angle \theta_a(r) vs Radius');
xlabel('Radius r (µm)');
ylabel('\theta_a (degrees)');

% Subplot 2: Refractive Index Profile
subplot(2,1,2);
hold on;
grid on;
title('b) Core Refractive Index n_1(r) vs Radius');
xlabel('Radius r (µm)');
ylabel('Refractive Index, n_1(r)');

for g = g_values
    % Core index profile
    n1_r = n1_0 * sqrt(1 - 2 * Delta * (r/a).^g);

    % Numerical Aperture
    NA_r = sqrt(max(n1_r.^2 - n2^2, 0)); % Avoid sqrt of negative values

    % Acceptance Angle in degrees
    theta_a = real(asind(min(NA_r, 1)));  % Clamp NA to 1 to avoid domain errors

    % --- CORRECTED PLOTTING ORDER ---
    
    % Plot theta_a(r) on the first subplot
    subplot(2,1,1);
    plot(r*1e6, theta_a, 'DisplayName', ['g = ' num2str(g)], 'LineWidth', 2);

    % Plot n1(r) on the second subplot
    subplot(2,1,2);
    plot(r*1e6, n1_r, 'DisplayName', ['g = ' num2str(g)], 'LineWidth', 2);
end

% Add legends
subplot(2,1,1);
legend('show');
subplot(2,1,2);
legend('show');