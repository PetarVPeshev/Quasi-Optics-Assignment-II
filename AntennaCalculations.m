close all;
clear;
clc;

%% Constants and Parameters
% Indecies
N = 300;                        % Number of sample points
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

%% Theta and Phi-Components of Spherical Cooridnates
th = linspace(eps, pi, N);
dth = th(2) - th(1);
ph = linspace(eps, 2 * pi, N);
dph = ph(2) - ph(1);
[ TH, PH ] = meshgrid(th, ph);

%% x, y, and z-Components of Wave Number
KX = k * sin(TH) .* cos(PH);
KY = k * sin(TH) .* sin(PH);
kz = -1j * sqrt( -(k^2 - KX.^2 - KY.^2) );

%% Calculate Spectral Green's Function (SGF)
ej_SGF = calculateEJ_SGF(er, k, KX, KY);
plotSGF(ej_SGF, k, KX, 'EJ');

%% Calculate Fourier Transform (FT) of Current Distribution
Jf = circFTCurrent(k, J0, TH, af, p);
plotAiryCurrent(Jf, TH);

%% Calculate Electric Far-Field of Feed
Ef = calculateEFarfield( ej_SGF, Jf, k, R, TH, PH, kz );
plotFarfield(Ef, TH, PH, N);
xlim([-0.25 0.25]);
ylim([-0.25 0.25]);
xticks(( -0.25 : 0.05 : 0.25 ));
yticks(( -0.25 : 0.05 : 0.25 ));
plotFarfield(Ef, TH, PH, N);
view(-130, 30);
xlim([-0.25 0.25]);
ylim([-0.25 0.25]);
xticks(( -0.25 : 0.05 : 0.25 ));
yticks(( -0.25 : 0.05 : 0.25 ));

%% Calculate Equivalent Aperture Current Distribution
[ J, M, RHO, PHc ] = calculatePACurrent( Ef, Z, k, R, TH, PH, rf2rD, rf );
plotCylCurrent(J, RHO, PHc);

%% Calculate Current Distribution Fourier Transform (FT)
tic
Jft = calculateCylFTCurrent( J, KX, KY, RHO, PHc );
toc

%% Calculate Electric Far-Field
Jft = convertCylToCar( Jft, PH );
E = calculateEFarfield( ej_SGF, Jft, k, R, TH, PH, kz );
plotFarfield(E, TH, PH, N);
caxis([-40, 0]);
view(-45, 20);
zlim([-40 0]);

%% Calculate Uniform Aperture Current Distribution Fourier Transform (FT)
Jun = circFTCurrent(k, J0, TH, rD / 2, p);
plotAiryCurrent(Jun, TH);
ylim([-100 0]);

%% Calculate Electric Far-Field of Uniform Aperture
Eun = calculateEFarfield( ej_SGF,  Jun, k, R, TH, PH, kz );
plotFarfield(Eun, TH, PH, N);
caxis([-40, 0]);
view(-25, 50);
zlim([-40 0]);

%% Calculate Antenna Efficiencies
[ Te, Se, Ae ] = calculateREfficiency( E, Ef, R, TH, PH, Z, k, rD, rf );

%% Calculate Maximum Possible Directivity, Directivity, and Gain
[ Dm, D, G ] = calculateRParameters( rD, wlen, Te, Ae );
