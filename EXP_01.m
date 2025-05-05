clc;
close all;
clear;

%Inputs..........
n1 = 1.48;      %Core Refractive Index
n2 = 0.46;      %Cladding Refractive Index
delta = 0.015;      %Relative Refractive index difference
lamda = 0.85e-6;       %Operating Wavelength
a = 40e-6;     %Radius of the core
k = 2*pi/lamda;     %Phase propagation constant
theta = 28.90;          %Acceptance angle

%Calculate Numerical Aperture(NA)
NA = n1*sqrt(2*delta);
disp('Numerical Aperture (NA): ')
disp(NA);

%Propagation Constant..............
beta = (2*pi*n1*cos(theta))/lamda;
disp('Propagation Constant: ');
disp(k);

%V Number................
V = (2*pi*a*NA)/lamda;
disp('V Number:' );
disp(V);

%Cutoff Wavelength..........
Vc = 2.405;     %Only for Single Mode Fiber
lamda_c = (2*pi*a*NA)/Vc;
disp('Cutoff Wavelength: ');
disp(lamda_c);

%(b)--------------------------------------------
%Check for single-mode or multi-mode
if V < 2.405
    disp('Single-mode fiber');
else
    disp('Multi-mode fiber');
end
disp(' ');

%Number of modes traveling in Fiber (Ms)
Ms = (V^2)/2;
disp('Number of modes traveling in Fiber (Ms):');
disp(Ms);

%Normalized Propagation Constant.............
% NPC = ((beta^2) - (n2^2)*(k^2))/(n1*(k^2) - (n2^2)*(k^2));
% disp('Normalized Propagation Constant: ');
% disp(NPC);

%Graphical representation of V verses normalized propogation
%constant................
lambda_range = linspace(1e-6, 1.8e-6, 100);  % Range of wavelengths

% Pre-allocate arrays for efficiency
b_values = zeros(size(lambda_range));
V_values = zeros(size(lambda_range));

%Calculate b and V for each wavelength
for i = 1:length(lambda_range)
    lambda = lambda_range(i);
    b_values(i) = sqrt((n1^2 - n2^2) * (2*pi/lambda)^2 * a^2);
    V_values(i) = 2*pi * a / lambda * sqrt(n1^2 - n2^2);
end

%Plot b vs. V
plot(V_values, b_values);
xlabel('V number (V)');
ylabel('Normalized propagation constant (b)');
title('Normalized Propagation Constant (b) vs. V Number (V)');
grid on;