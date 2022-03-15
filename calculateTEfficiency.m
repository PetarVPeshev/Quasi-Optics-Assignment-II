function [ Te ] = calculateTEfficiency( E, k, R, TH, PH, RHO, D, f )
%calculateTEfficiency This function calculates the taper efficiency of a
%reflector antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Element Step
    dph = PH(2, 1) - PH(1, 1);
    drho = RHO(1, 2) - RHO(1, 1);
    %% Calculate Aperture Field
    E = E * 2 * pi * R / exp(-1j * k * R);
    Ea = zeros( size(E) );
    Ea(:, :, 1) = - E(:, :, 3) .* cos(PH) .* ( cos( TH / 2 ) .^ 2) * ...
                  exp(-1j * k * 2 * f) / f;
    Ea(:, :, 2) = - E(:, :, 2) .* ( cos( TH / 2 ) .^ 2) * ...
                  exp(-1j * k * 2 * f) / f;
    %% Calculate Aperture Far-Field Magnitude
    Eat = abs( Ea(:, :, 1) ).^2 + abs( Ea(:, :, 2) ).^2 + ...
         abs( Ea(:, :, 3) ).^2;
    %% Calculate Effective Antenna Area
    Exi = sum( sum( abs( Ea(:, :, 1) ) .* RHO ) ) * drho * dph;
    Eyi = sum( sum( abs( Ea(:, :, 2) ) .* RHO ) ) * drho * dph;
    Ezi = sum( sum( abs( Ea(:, :, 3) ) .* RHO ) ) * drho * dph;
    Eam = ( Exi + Eyi + Ezi ) ^ 2;
    Eai = sum( sum( Eat .* RHO ) ) * drho * dph;
    Aeff = Eam / Eai;
    %% Calculate Taper Efficiency
    A = pi * (D^2) / 4;
    Te = Aeff / A;
    Te( Te > 1 ) = 1;
    Te = Te * 100;
end