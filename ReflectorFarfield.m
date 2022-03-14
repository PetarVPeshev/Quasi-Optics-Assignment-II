close all;
clear;
clc;

%% Constants and Parameters
% Indecies
N = 100;                        % Number of sample points
% Field Parameters
f = 100 * 1e9;                  % Source frequency [Hz]
R = 1;                          % Radial distance [m]
% Antenna Parameters
Df2wlen = 4;                    % Feed diameter to wavelength
J0 = 1;                         % Amplitude of feed current distribution
p = [0 1 0];                    % Polarization of feed current distribution
rf = 1;                         % Antenna focal length [m]
rf2rD = 5;                      % Antenna focal length to diameter
% Medium
er = 1;                         % Relative permittivity
c = physconst('LightSpeed');    % Speed of light [m/s]
e0 = 8.8541878128 * 1e-12;      % Permittivity of free space [F/m]
u0 = 4 * pi * 1e-7;             % Permeability of free space [H/m]

%% Parameters
wlen = c / f;                   % Wavelength [m]
k = 2*pi / wlen;                % Magnitude of wave number [rad/m]
Z = sqrt( u0 / (e0 * er) );     % Wave impedance [Ohm]
af = Df2wlen * wlen / 2;        % Radius of feed [m]
rD = rf / rf2rD;                % Diameter of reflector [m]

%% Rho and Phi-Components of Cylindrical Coordinates
ph = linspace(eps, 2 * pi, N);
dph = ph(2) - ph(1);
rho = linspace(eps, rD / 2, N);
drho = rho(2) - rho(1);
[ RHO, PH ] = meshgrid(rho, ph);
TH = 2 * atan( RHO / (2 * rf) );

%% x, y, and z-Components of Wave Number
KX = k * sin(TH) .* cos(PH);
KY = k * sin(TH) .* sin(PH);
kz = -1j * sqrt( -(k^2 - KX.^2 - KY.^2) );

%% Calculate FT of Feed Current Distribution
Jf = circFTCurrent( k, J0, TH, af, p );

%% Calculate New Spectral Green's Function (SGF)
ej_SGF = calculateEJ_SGF( er, k, KX, KY );

%% Calculate Electric Far-Field of Feed
Ef = calculateEFarfield( ej_SGF, Jf, k, R, TH, kz );
Ef( isnan(Ef) ) = 0;

%% Calculate Equivalent Aperture Current Distribution
[ J, M ] = calculatePACurrent( Ef, Z, k, R, TH, rf );

%% Calculate Current Distribution Fourier Transform (FT)
Jft = calculateCylFTCurrent( J, KX, KY, RHO, PH );

%% Calculate Electric Far-Field
E = calculateEFarfield( ej_SGF, Jft, k, R, TH, kz );
E( isnan(E) ) = 0;
plotCarFarfield(E, RHO, PH);
caxis([-40, 0]);
zlim([-150 0]);

%% Calculate Uniform Aperture Current Distribution Fourier Transform (FT)
Jun = circFTCurrent(k, J0, TH, rD / 2, p);

%% Calculate Electric Far-Field of Uniform Aperture
Eun = calculateEFarfield( ej_SGF, Jun, k, R, TH, kz );
Eun( isnan(Eun) ) = 0;
plotCarFarfield(Eun, RHO, PH);
caxis([-40, 0]);
zlim([-150 0]);

%% Plot in 1D
% Define Theta from - Theta_max to Theta_max
th = zeros( 1, 2 * size(TH, 2) );
th( size(TH, 2) + 1 : end ) = TH(1, :);
th( 1 : size(TH, 2) ) = - rot90( TH(1, :), 2 );
% Extract E field magnitude and uniform E field
Eth = zeros( 1, size(th, 2) );
Eth( size(E, 2) + 1 : end ) = sqrt( abs( E(51, :, 1) ).^2 + ...
                           abs( E(51, :, 2) ).^2 + abs( E(51, :, 3) ).^2 );
Eth( 1 : size(E, 2) ) = rot90(sqrt( abs( E(1, :, 1) ).^2 + ...
                     abs( E(1, :, 2) ).^2 + abs( E(1, :, 3) ).^2 ), 2);
Euth = zeros( 1, size(th, 2) );
Euth( size(E, 2) + 1 : end ) = sqrt( abs( Eun(51, :, 1) ).^2 + ...
                       abs( Eun(51, :, 2) ).^2 + abs( Eun(51, :, 3) ).^2 );
Euth( 1 : size(E, 2) ) = rot90(sqrt( abs( Eun(1, :, 1) ).^2 + ...
                     abs( Eun(1, :, 2) ).^2 + abs( Eun(1, :, 3) ).^2 ), 2);
% Plot
figure();
plot(th * 180 / pi, 20 * log10( Eth ) - max( 20 * log10( Euth ) ), ...
     'LineWidth', 3.0);
hold on;
plot(th * 180 / pi, 20 * log10( Euth ) - max( 20 * log10( Euth ) ), ...
     '--', 'LineWidth', 3.0);
grid on;
xlabel('\theta [deg]');
ylabel('|E| [dB] at XZ Plane');
xlim([min(th * 180 / pi) max(th * 180 / pi)]);
ylim([-60 0]);
legend('|E|', '|E_{JU}|');
