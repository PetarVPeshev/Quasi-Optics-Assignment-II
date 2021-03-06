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

%% x, y, and z-Components of Wave Number for Feed
KXf = k * sin(TH) .* cos(PH);
KYf = k * sin(TH) .* sin(PH);
kzf = -1j * sqrt( -(k^2 - KXf.^2 - KYf.^2) );

%% Calculate Spectral Green's Function (SGF) for Feed
ej_SGF_f = calculateEJ_SGF( er, k, KXf, KYf );
plotSGF(ej_SGF_f, k, KXf, 'EJ');

%% Calculate Fourier Transform (FT) of Current Distribution
Jf = circFTCurrent( k, J0, TH, af, p );
plotAiryCurrent(Jf, TH);

%% Calculate Electric Far-Field of Feed
Ef = calculateEFarfield( ej_SGF_f, Jf, k, R, TH, kzf );

%% Convert to Spherical Coordinates
Ef = convertCarToSph( Ef, TH, PH );
plotFarfield( Ef, TH, PH );
xlim([-0.25 0.25]);
ylim([-0.25 0.25]);
xticks(( -0.25 : 0.05 : 0.25 ));
yticks(( -0.25 : 0.05 : 0.25 ));
plotFarfield( Ef, TH, PH );
view(-130, 30);
xlim([-0.25 0.25]);
ylim([-0.25 0.25]);
xticks(( -0.25 : 0.05 : 0.25 ));
yticks(( -0.25 : 0.05 : 0.25 ));

%% Rho and Phi-Components of Cylindrical Coordinates
rho = linspace(eps, rD / 2, N);
[ RHO, PHP ] = meshgrid(rho, ph);
THP = 2 * atan( RHO / (2 * rf) );
rhoi = eps : ( rho(2) - rho(1) ) : ( 2 * rf * tan( pi / 4 ) );
[ RHOi, PHPi ] = meshgrid(rhoi, ph);
THPi = 2 * atan( RHOi / (2 * rf) );

%% Calculate FT of Feed Current Distribution in New Coordinates
Jff = circFTCurrent( k, J0, THP, af, p );

%% x, y, and z-Components of Wave Number in New Coordinates
KX = k * sin(THP) .* cos(PHP);
KY = k * sin(THP) .* sin(PHP);
kz = -1j * sqrt( -(k^2 - KX.^2 - KY.^2) );
KXi = k * sin(THPi) .* cos(PHPi);
KYi = k * sin(THPi) .* sin(PHPi);
kzi = -1j * sqrt( -(k^2 - KXi.^2 - KYi.^2) );

%% Calculate New Spectral Green's Function (SGF) in New Coordinates
ej_SGF = calculateEJ_SGF( er, k, KX, KY );
plotSGF( ej_SGF, k, KX, 'EJ' );

%% Calculate Electric Far-Field of Feed in New Coordinates
Eff = calculateEFarfield( ej_SGF, Jff, k, R, THP, kz );

%% Calculate Equivalent Aperture Current Distribution
[ J, M ] = calculatePACurrent( Eff, Z, k, R, THP, PHP, rf );
plotCurrent(M, RHO, PHP, 'M');
plotCurrent(J, RHO, PHP, 'J');

%% Calculate Current Distribution Fourier Transform (FT)
tic
Jft = calculateCylFTCurrent( J, KX, KY, RHO, PHP );
toc

%% Calculate Electric Far-Field
E = calculateEFarfield( ej_SGF, Jft, k, R, THP, kz );
plotCarFarfield(E, RHO, PHP);
caxis([-40, 0]);
view(-45, 20);
zlim([-40 0]);

%% Calculate Uniform Aperture Current Distribution Fourier Transform (FT)
Jun = circFTCurrent(k, J0, THP, rD / 2, p);

%% Calculate Electric Far-Field of Uniform Aperture
Eun = calculateEFarfield( ej_SGF, Jun, k, R, THP, kz );
plotCarFarfield(Eun, RHO, PHP);
caxis([-40, 0]);
view(-25, 50);
zlim([-40 0]);

%% Calculate Electric Far-Field of Feed in 0 to 90 Degrees
Jfi = circFTCurrent( k, J0, THPi, af, p );
ej_SGFi = calculateEJ_SGF( er, k, KXi, KYi );
Efi = calculateEFarfield( ej_SGFi, Jfi, k, R, THPi, kzi );

%% Calculate Antenna Efficiencies
Se = calculateSOEfficiency(Eff, Efi, Z, k, R, THP, PHP, THPi, PHPi, rf);
Te = calculateTEfficiency(Eff, k, R, THP, PHP, RHO, rD, rf);
Se = Se / 100;
Te = Te / 100;
Ae = Se * Te;

%% Calculate Maximum Possible Directivity, Directivity, and Gain
[ Dm, D, G ] = calculateRParameters( rD, wlen, Te, Ae );
