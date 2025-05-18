%% Hohmann Transfer from Earth to Venus (3-Dimension)
clear;
clc;

%% Constants
mu_Sun = 1.32712440018e11;  % Gravitational parameter of the Sun (km^3/s^2)
AU = 1.496e8;                % Astronomical unit (km)
shift_distance = 0.1e8;      % Shift distance along z-axis (km)

%% Orbits
r_earth = 1*AU;             % Radius of Earth's orbit (km)
r_venus = 0.72*AU;          % Radius of Venus' orbit (km)

%% Hohmann Transfer Calculations

% Semi-major axis of transfer ellipse
a_transfer = (r_earth + r_venus) / 2;

% Eccentricity of the transfer ellipse
e = (r_venus - r_earth) / (r_venus + r_earth);

% Velocity at Earth's orbit
v_earth = sqrt(mu_Sun / r_earth);

% Velocity at transfer orbit (at the point of departure from Earth)
v_transfer_departure = sqrt((2 * mu_Sun / r_earth) - (mu_Sun / a_transfer));

% Velocity at Venus' orbit
v_venus = sqrt(mu_Sun / r_venus);

% Velocity at transfer orbit (at the point of arrival at Venus)
v_transfer_arrival = sqrt((2 * mu_Sun / r_venus) - (mu_Sun / a_transfer));

% Delta-v required for transfer
delta_v_departure = v_transfer_departure - v_earth;
delta_v_arrival = v_venus - v_transfer_arrival;

%% Plot Trajectories
theta = linspace(0, 2*pi, 1000);

% Earth's orbit
xE = r_earth * cos(theta);
yE = r_earth * sin(theta);
zE = zeros(size(xE))+ 2e8;

% Venus' orbit
xV = r_venus * cos(theta);
yV = r_venus * sin(theta);
zV = zeros(size(xV)) + 2e8;

% Transfer orbit
theta1 = linspace(0, pi, 1000);
r_transfer = (a_transfer * (1 - e^2)) ./ (1 + e * cos(theta1));
xT = r_transfer .* cos(theta1);
yT = r_transfer .* sin(theta1);
zT = zeros(size(xT)) + 2e8;

figure;
plot3(xE, yE, zE, 'b'); hold on;
plot3(xV, yV, zV, 'g');
plot3(xT, yT, zT, 'r--');
plot3(0, 0, 2e8, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', 'y'); % Sun
axis equal;
legend('Earth Orbit', 'Venus Orbit', 'Transfer Orbit', 'Sun');
title('Hohmann Transfer from Earth to Venus');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');

% Add some perspective to the z-axis
zlim([0, 4e8]);
zticks(0:0.5e8:4e8);
%% Display Delta-v
disp(['Delta-v required for departure from Earth: ' num2str(delta_v_departure) ' km/s']);
disp(['Delta-v required for arrival at Venus: ' num2str(delta_v_arrival) ' km/s']);

