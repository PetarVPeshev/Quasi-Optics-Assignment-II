close all;
clear;
clc;

%% Constants and Parameters
% Indecies
N = 500;                        % Number of sample points
% Field Parameters
f = 100 * 1e9;                  % Source frequency [Hz]
R = 1;                          % Radial distance [m]
% Antenna Parameters
Df2wlen = 4;                    % Feed diameter to wavelength
J0 = 1;                         % Amplitude of feed current distribution
p = [0 1 0];                    % Polarization of feed current distribution
rf = 1;                         % Antenna focal length [m]
rf2rD = 0.5;                    % Antenna focal length to diameter
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

%% Calculate Equivalent Aperture Current Distribution
[ J, M ] = calculatePACurrent( Ef, Z, k, R, TH, PH, rf );
plotCurrent(M, RHO, PH, 'M');
