function [ Te, Se, Ae ] = calculateREfficiency( E, Ef, Efi, R, TH, PH, RHO, THi, PHi, Z, k, D, f )
%calculateREfficiency This function calculates the taper, spill-over, and
%aperture efficiency of parabolic antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Step Sizes
    dth = TH(1, 2) - TH(1, 1);
    dph = PH(2, 1) - PH(1, 1);
    drho = RHO(1, 2) - RHO(1, 1);
    %% Calculate Far-Field Magnitude
    % Electric far-field of feed in 0 to Theta max region
    Ef = Ef * R * exp(-1j * k * f) / ( f * exp(-1j * k * R) );
    Efsph = convertCarToSph( Ef, TH, PH );
    Eft = abs( Efsph(:, :, 1) ).^2 + abs( Efsph(:, :, 2) ).^2 + ...
          abs( Efsph(:, :, 3) ).^2;
    % Electric far-field of feed in 0 to 90 Degrees region
    Efi = Efi * R * exp(-1j * k * f) / (f * exp(-1j * k * R));
    Efisph = convertCarToSph( Efi, THi, PHi );
    Efti = abs( Efisph(:, :, 1) ).^2 + abs( Efisph(:, :, 2) ).^2 + ...
           abs( Efisph(:, :, 3) ).^2;
    % Electric far-field of antenna
    Et = abs( E(:, :, 1) ).^2 + abs( E(:, :, 2) ).^2 + ...
         abs( E(:, :, 3) ).^2;
    %% Calculate Radiation Intensity of Feed Far-Field
    Uf = Eft * 4 * (f .^ 2) / (2 * Z);
    Ufi = Efti * 4 * (f .^ 2) / (2 * Z);
    %% Remove Singularities
    Ufi( isnan(Ufi) ) = 0;
    %% Calculate Radiated Power of Feed
    Pf = sum( sum( Uf .* sin(TH) ) ) * dth * dph;
    Pfi = sum( sum( Ufi .* sin(THi) ) ) * dth * dph;
    %% Calculate Spill-Over Efficiency
    Se = Pf / Pfi;
    %% Calculate Effective Antenna Area
    Eam = ( ( sum( sum( abs( E(:, :, 1) ) .* RHO ) ) * drho * dph ) + ...
               ( sum( sum( abs( E(:, :, 2) ) .* RHO ) ) * drho * dph ) + ...
               ( sum( sum( abs( E(:, :, 3) ) .* RHO ) ) * drho * dph ) ) .^ 2;
    Ea = sum( sum( Et .* RHO ) ) * drho * dph;
    Aeff = Eam / Ea;
    %% Calculate Taper Efficiency
    A = pi * (D / 2)^2;
    Te = Aeff / A;
    %% Calculate Aperture Efficiency and Convert to Percentage
    Ae = Se * Te * 100;
    Se = Se * 100;
    Te = Te * 100;
end