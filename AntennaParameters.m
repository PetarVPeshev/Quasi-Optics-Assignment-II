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
rf2rD = linspace(0.5, 5, 20);   % Antenna focal length to diameter
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
rD = rf ./ rf2rD;            % Diameter of reflector [m]

%% Phi-Components of Cylindrical Coordinates
ph = linspace(eps, 2 * pi, N);
dph = ph(2) - ph(1);

%% Calculate Antenna Parameters
Te = zeros(1, size(rf2rD, 2) );
Se = zeros(1, size(rf2rD, 2) );
Ae = zeros(1, size(rf2rD, 2) );
Dm = zeros(1, size(rf2rD, 2) );
D = zeros(1, size(rf2rD, 2) );
G = zeros(1, size(rf2rD, 2) );
for i = 1:size(rf2rD, 2)
    %% Rho and Theta-Components of Cylindrical Coordinates
    rho = linspace(eps, rD(i) / 2, N);
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
    
    %% Rho and Phi-Components of Cylindrical Coordinates in 0 to 90 Degrees
    rhoi = eps : drho : ( 2 * rf * tan( pi / 4 ) );
    [ RHOi, PHi ] = meshgrid(rhoi, ph);
    THi = 2 * atan( RHOi / (2 * rf) );
    
    %% x, y, and z-Components of Wave Number in 0 to 90 Degrees
    KXi = k * sin(THi) .* cos(PHi);
    KYi = k * sin(THi) .* sin(PHi);
    kzi = -1j * sqrt( -(k^2 - KXi.^2 - KYi.^2) );
    
    %% Calculate Electric Far-Field of Feed in 0 to 90 Degrees
    Jfi = circFTCurrent( k, J0, THi, af, p );
    ej_SGFi = calculateEJ_SGF( er, k, KXi, KYi );
    Efi = calculateEFarfield( ej_SGFi, Jfi, k, R, THi, kzi );
    
    %% Calculate Antenna Efficiencies
    [ Te(i), Se(i), Ae(i) ] = calculateREfficiency( E, Ef, Efi, R, TH, PH, RHO, ...
                                           THi, Z, k, rD(i), rf );
    
    %% Calculate Maximum Possible Directivity, Directivity, and Gain
    [ Dm(i), D(i), G(i) ] = calculateRParameters( rD(i), wlen, Te(i), Ae(i) );
end

%% Plots
% Efficiencies
figure();
plot(rf2rD, Te, 'LineWidth', 3.0);
hold on;
plot(rf2rD, Se, 'LineWidth', 3.0);
hold on;
plot(rf2rD, Ae, 'LineWidth', 3.0);
grid on;
xlabel('Focal length to diameter ratio [f / D]');
ylabel('%');
legend('\eta_{T}', '\eta_{S}', '\eta_{A}');
% Parameters
figure();
plot(rf2rD, 20 * log10(Dm), 'LineWidth', 3.0);
hold on;
plot(rf2rD, 20 * log10(D), 'LineWidth', 3.0);
hold on;
plot(rf2rD, 20 * log10(G), 'LineWidth', 3.0);
grid on;
xlabel('Focal length to diameter ratio [f / D]');
ylabel('[dB]');
legend('D_{max}','D','G');
