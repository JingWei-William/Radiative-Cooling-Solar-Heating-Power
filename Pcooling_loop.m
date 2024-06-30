% Constants
h = 6.626e-34; % Planck's constant (J·s)
c = 3e8; % Speed of light (m/s)
kB = 1.381e-23; % Boltzmann constant (J/K)

% Parameters
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

% Interpolate solar spectrum data to match lambda range
IAM1_5_interp = interp1(lambda_sun, IAM1_5, lambda, 'linear', 0);

% Blackbody spectral radiance function
IBB = @(T, lambda) (2 * h * c^2 ./ lambda.^5) ./ (exp(h * c ./ (lambda * kB * T)) - 1);

% Power radiated out by the cooler function
Prad_func = @(T) trapz(theta, 2 * pi * sin(theta) .* cos(theta)) * trapz(lambda, IBB(T, lambda) .* epsilon_surf_interp);

% Incident atmospheric thermal radiation function
Patm_func = trapz(theta, 2 * pi * sin(theta) .* cos(theta)) * trapz(lambda, IBB(Tamb, lambda) .* epsilon_atm_interp);

% Incident solar power absorbed by the structure function
PSun_func = trapz(lambda, IAM1_5_interp .* epsilon_surf_interp);

% Define the range of temperatures to calculate
T_values = 288:1:338; % Example range from 270K to 330K in steps of 10K

% Initialize results storage
Prad_values = zeros(size(T_values));
Patm_values = zeros(size(T_values));
PSun_values = zeros(size(T_values));
Pcond_conv_values = zeros(size(T_values));
Pcool_density_values = zeros(size(T_values));
Pcool_values = zeros(size(T_values));

% Calculate results for each temperature
for i = 1:length(T_values)
    T = T_values(i);
    
    % Power radiated out by the cooler
    Prad = Prad_func(T);
    
    % Incident atmospheric thermal radiation
    Patm = Patm_func;
    
    % Incident solar power absorbed by the structure
    PSun = PSun_func;
    
    % Power lost due to convection and conduction
    Pcond_conv = hc * (Tamb - T);
    
    % Net cooling power density
    Pcool_density = Prad - Patm - PSun - Pcond_conv;
    
    % Total net cooling power
    Pcool = Pcool_density * A;
    
    % Store results in arrays
    Prad_values(i) = Prad;
    Patm_values(i) = Patm;
    PSun_values(i) = PSun;
    Pcond_conv_values(i) = Pcond_conv;
    Pcool_density_values(i) = Pcool_density;
    Pcool_values(i) = Pcool;
end

% Display results
fprintf('%6s %10s %10s %10s %15s %10s\n', 'Temp(K)', 'Prad(W/m^2)', 'Patm(W/m^2)', 'PSun(W/m^2)', 'Pcool_density(W/m^2)', 'Pcool(W)');
for i = 1:length(T_values)
    fprintf('%6.0f %10.2f %10.2f %10.2f %15.2f %10.2f\n', T_values(i), Prad_values(i), Patm_values(i), PSun_values(i), Pcool_density_values(i), Pcool_values(i));
end
