function [ Te, Se, Ae ] = calculateREfficiency( E, Ef, R, TH, PH, Z, k, D, f )
%calculateREfficiency This function calculates the taper, spill-over, and
%aperture efficiency of parabolic antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Maximum Elevation Angle
    fsh = f / D;
    th0 = 2 * acot( 4 * fsh );
    %% Calculate Elevation Angles on Reflector
    THe = TH - pi / 2;
    THe( abs( TH - pi / 2 ) > th0 ) = 0;
    %% Convert Electric Field to Cartesian Coordinates
    RHO = R .* sin(TH);
    drho = RHO(1, 2) - RHO(1, 1);
    dph = PH(2, 1) - PH(1, 1);
    E = convertSphToCyl(E, TH);
    %% Calculate Far-Field Magnitude
    Ef = Ef * R * exp(-1j * k * f) / ( f * exp(-1j * k * R) );
    Eft = abs( Ef(:, :, 1) ).^2 + abs( Ef(:, :, 2) ).^2 + ...
          abs( Ef(:, :, 3) ).^2;
    Et = abs( E(:, :, 1) ).^2 + abs( E(:, :, 2) ).^2 + ...
         abs( E(:, :, 3) ).^2;
    %% Calculate Radiation Intensity of Feed Far-Field
    Uf = Eft * (f .^ 2) / (2 * Z);
    %% Calculate Spillover Efficiency
    Se = sum( sum( Uf .* sin(THe) ) ) / ...
         sum( sum( Uf .* sin(TH - pi / 2) ) );
    %% Calculate Taper Efficiency
    Aeff = ( abs( sum( sum( E(:, :, 1) .* RHO ) ) * drho * dph + ...
           sum( sum( E(:, :, 2) .* RHO ) ) * drho * dph + ...
           sum( sum( E(:, :, 3) .* RHO ) ) * drho * dph ) .^2 ) / ...
           ( sum( sum( Et .* RHO ) ) * drho * dph );
    A = pi * (D / 2)^2;
    Te = Aeff / A;
    %% Calculate Aperture Efficiency
    Ae = Se * Te;
end