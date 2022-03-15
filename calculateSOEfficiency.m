function [ Se ] = calculateSOEfficiency( E, Ei, Z, k, R, TH, PH, THi, PHi, f)
%calculateSOEfficiency This function calculates the spill-over efficiency
%of a reflector antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Element Step
    dth = TH(1, 2) - TH(1, 1);
    dph = PH(2, 1) - PH(1, 1);
    %% Calculate Far-Field Magnitude
    % Electric far-field of feed in 0 to Theta max region
    E = E * R * exp(-1j * k * f) / ( f * exp(-1j * k * R) );
    Esph = convertCarToSph( E, TH, PH );
    Et = abs( Esph(:, :, 1) ).^2 + abs( Esph(:, :, 2) ).^2 + ...
          abs( Esph(:, :, 3) ).^2;
    % Electric far-field of feed in 0 to 90 Degrees region
    Ei = Ei * R * exp(-1j * k * f) / (f * exp(-1j * k * R));
    Eisph = convertCarToSph( Ei, THi, PHi );
    Eit = abs( Eisph(:, :, 1) ).^2 + abs( Eisph(:, :, 2) ).^2 + ...
           abs( Eisph(:, :, 3) ).^2;
    %% Calculate Radiation Intensity of Feed Far-Field
    Uf = Et * 4 * (f .^ 2) / (2 * Z);
    Ufi = Eit * 4 * (f .^ 2) / (2 * Z);
    %% Remove Singularities
    Ufi( isnan(Ufi) ) = 0;
    %% Calculate Radiated Power of Feed
    Pf = sum( sum( Uf .* sin(TH) ) ) * dth * dph;
    Pfi = sum( sum( Ufi .* sin(THi) ) ) * dth * dph;
    %% Calculate Spill-Over Efficiency
    Se = Pf / Pfi;
    Se = Se * 100;
end