% Constants
h = 6.626e-34; % Planck's constant (J·s)
c = 3e8; % Speed of light (m/s)
kB = 1.381e-23; % Boltzmann constant (J/K)

% Parameters
T = 288; % Temperature of the radiative cooler (K)
Tamb = 298; % Ambient temperature (K)
hc = 6.9; % Non-radiative heat transfer coefficient (W/m^2·K)
A = 0.0009; % Area of the radiative cooler (m^2)

% Load solar spectrum data
solar_spectrum = xlsread('Solar Radiation Spectrum.xlsx');
lambda_sun = solar_spectrum(:, 1); % Wavelength in micrometers (um)
IAM1_5 = solar_spectrum(:, 2); % Solar spectrum (W/m^2/um)

% Convert wavelength from um to m
lambda_sun = lambda_sun * 1e-6; % Convert wavelength from um to m

% Convert solar spectrum from W/m^2/um to W/m^2/m
IAM1_5 = IAM1_5 * 1e6; % Convert solar spectrum units

% Load atmospheric emissivity data
atmospheric_emissivity = xlsread('Atmospheric Window.xlsx');
lambda_atm = atmospheric_emissivity(:, 1) * 1e-6; % Convert wavelength from um to m
epsilon_atm_data = atmospheric_emissivity(:, 2); % Atmospheric emissivity

% Load surface emissivity data
surface_emissivity = xlsread('RCC Surface Emissivity.xlsx');
lambda_surf = surface_emissivity(:, 1) * 1e-6; % Convert wavelength from um to m
epsilon_surf_data = surface_emissivity(:, 2); % Surface emissivity

% Spectral and Angular Emissivity (example values, need to be provided)
lambda = linspace(0.3e-6, 20e-6, 1000); % Wavelength range (m)
theta = linspace(0, pi/2, 100); % Angle range (radians)

% Interpolate atmospheric emissivity data to match lambda range
epsilon_atm_interp = interp1(lambda_atm, epsilon_atm_data, lambda, 'linear', 0);

% Interpolate surface emissivity data to match lambda range
epsilon_surf_interp = interp1(lambda_surf, epsilon_surf_data, lambda, 'linear', 0);

% Blackbody spectral radiance
IBB = @(T, lambda) (2 * h * c^2 ./ lambda.^5) ./ (exp(h * c ./ (lambda * kB * T)) - 1);

% Power radiated out by the cooler
dV = 2 * pi * sin(theta) .* cos(theta);
Prad = trapz(theta, dV) * trapz(lambda, IBB(T, lambda) .* epsilon_surf_interp);

% Incident atmospheric thermal radiation
Patm = trapz(theta, dV) * trapz(lambda, IBB(Tamb, lambda) .* epsilon_atm_interp);

% Interpolate solar spectrum data to match lambda range
IAM1_5_interp = interp1(lambda_sun, IAM1_5, lambda, 'linear', 0);

% Incident solar power absorbed by the structure
PSun = trapz(lambda, IAM1_5_interp .* epsilon_surf_interp);

% Power lost due to convection and conduction
Pcond_conv = hc * (Tamb - T);

% Net cooling power density
Pcool_density = Prad - Patm - PSun - Pcond_conv;

% Total net cooling power
Pcool = Pcool_density * A;

% Display results
fprintf('Power radiated out by the cooler: %.2f W/m^2\n', Prad);
fprintf('Incident atmospheric thermal radiation: %.2f W/m^2\n', Patm);
fprintf('Incident solar power absorbed by the structure: %.2f W/m^2\n', PSun);
fprintf('Power lost due to convection and conduction: %.2f W/m^2\n', Pcond_conv);
fprintf('Net cooling power density: %.2f W/m^2\n', Pcool_density);
fprintf('Total net cooling power: %.2f W\n', Pcool);
