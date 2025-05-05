clc;
clear;
close all;

% Given parameters
a = 40e-6;      % Core radius in meters
n1 = 1.48;      % Core refractive index
n2 = 1.46;      % Cladding refractive index
lamda = 1.55e-6;        % Operating wavelength in meters

%Numerical Aperture (NA)
NA = sqrt(n1^2 - n2^2);

%V Number
V = 2 * pi * a * NA / lamda;

%Cutoff wavelength
Vc = 10.50;     %Considering multimode fiber........
lambda_c = (2*pi*a*NA)/Vc;

%Check whether the fiber is single mode or multimode
if V < 2.405
    fiber_type = 'Single mode';
else
    fiber_type = 'Multimode';
end

%Number of modes traveling in fiber (Ms)
Ms = (V^2) / 4;

%Graph between V number and wavelength
wavelength_range = linspace(lamda * 0.8, lamda * 1.2, 100);         % Varying wavelength
V_numbers = 2 * pi * a * NA ./ wavelength_range;
plot(wavelength_range * 1e6, V_numbers);
xlabel('Wavelength (\mum)');
ylabel('V Number');
title('V Number vs. Wavelength');
grid on;

%Display results
fprintf('Numerical Aperture (NA): %.4f\n', NA);
fprintf('V Number: %.4f\n', V);
fprintf('Cutoff wavelength: %.4f um\n', lambda_c* 1e6);
fprintf('Fiber type: %s\n', fiber_type);
fprintf('Number of modes traveling in fiber (Ms): %.4f\n', Ms);