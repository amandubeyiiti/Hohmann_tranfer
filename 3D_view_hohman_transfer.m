%% Hohmann Transfer from Earth to Venus (2-Dimension)
clear;
clc;
%% Define constants
mu_Sun = 1.32712440018e11;  
% Astronomical unit (km)
AU = 1.496e8;    

%% Orbital Parameters (km)
a_Earth = 1*AU; % radius of Earth's orbit around the Sun
a_venus = 0.72*AU; % radius of venus orbit around the Sun
%% Hohmann Transfer Calculations

% Velocities in Circular Orbits
v_Earth = sqrt(mu_Sun / a_Earth); % Current velocity in Earth orbit
v_venus = sqrt(mu_Sun / a_venus); % Velocity in Venus orbit


% Velocities in Transfer Orbit
a_transfer = (a_Earth + a_venus) / 2;
v_transfer_1 = sqrt(2 * mu_Sun / a_Earth - mu_Sun / a_transfer); % At departure from Earth
v_transfer_2 = sqrt(2 * mu_Sun / a_venus - mu_Sun / a_transfer); % At arrival at venus


% Delta-v Calculations
delta_v1 = abs(v_transfer_1 - v_Earth); % Needed to enter Hohmann transfer
delta_v2 = abs(v_transfer_2 - v_venus); % Needed to arrival at venus

% Calculate transfer orbit parameters
e_transfer = (a_Earth - a_venus) / (a_venus + a_Earth); % Eccentricity of transfer orbit
b_transfer = a_transfer * sqrt(1 - e_transfer^2); % Semi-minor axis

% Output Results
fprintf('Hohmann Transfer from Earth to venus:\n');
fprintf('Semi-major axis of transfer orbit: %.0f km\n', a_transfer)
fprintf('Delta-v for first burn at Earth: %.2f km/s\n', delta_v1)
fprintf('Delta-v for second burn at venus: %.2f km/s\n', delta_v2)
fprintf('Total Delta-v required: %.2f km/s\n', delta_v1 + delta_v2)

% For plotting of orbits

% Angles for plotting orbits
theta = linspace(0, 2 * pi, 1000);

% Earth Orbit (circle for simplification)
xEarth = a_Earth * cos(theta);
yEarth = a_Earth * sin(theta);

% Venus Orbit (circle for simplification)
xVenus = a_venus * cos(theta);
yVenus = a_venus * sin(theta);

% Transfer Orbit (ellipse)
theta_transfer = linspace(0, pi, 1000); 
xTransfer = a_transfer * cos(theta_transfer) - (a_transfer * e_transfer); 
yTransfer = b_transfer * sin(theta_transfer);

% Create figure for plotting
figure;
hold on;
grid on;
axis equal;
% Plot orbits
plot(xEarth, yEarth, 'b', 'LineWidth', 1); % Earth orbit in blue
plot(xVenus, yVenus, 'g', 'LineWidth', 1); % Venus orbit in red
plot(xTransfer, yTransfer, 'r--', 'LineWidth', 1.5); % Transfer orbit in dashed black
plot(0, 0, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', 'y'); % Sun
% Annotations and labels
title('Hohmann Transfer from Earth to Venus');
xlabel('Distance from Sun (km)');
ylabel('Distance from Sun (km)');
legend('Earth Orbit', 'Venus Orbit', 'Transfer Orbit', "sun");

