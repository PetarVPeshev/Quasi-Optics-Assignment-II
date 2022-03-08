function [ J, M, rho, PH ] = calculatePACurrent( Ef, Z, k, R, TH, PH, fsh, f )
%calculatePACurrent This function calculates the equivalent aperture
%current of parabolic antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
%   Note: optimize neglection of far-field out of reflector boundaries
    %% Calculate Maximum Elevation Angle
    th0 = 2 * acot( 4 * fsh );
    %% Calculate Valid Elevation Angles and Resize Polar Angles
%     THe = TH - pi / 2;
    TH = TH - pi / 2;
    TH( abs( TH ) > th0 ) = 0;
    THe = TH( : , any(TH) );
    PH = PH( :, any(TH) );
    %% Calculate Polar Distance
    rho = 2 * f * tan( THe / 2 );
    %% Calculate Electric Field Components with no Spherical Wave Part
    Ef = Ef * 2 * pi * R / exp(-1j * k * R);
    %% Neglect Electric Far-Field out of Reflector Boundaries
%     E = Ef;
    E = zeros([ size(THe), 3 ]);
    for i = 1:3
        temp = Ef(:, :, i);
        temp( abs( TH - pi / 2 ) > th0 ) = 0 + 0j;
        temp = temp( :, any(TH) );
        E(:, :, i) = temp;
    end
    %% Calculate Magnetic Current Distribution in Cylindrical Coordinates
    M = zeros( size(E) );
    M(:, :, 1) = - E(:, :, 3) .* ( cos( THe / 2 ) .^2 ) * ...
                 exp(-1j * k * 2 * f) / f;
    M(:, :, 2) = E(:, :, 2) .* ( cos( THe / 2 ) .^2 ) * ...
                 exp(-1j * k * 2 * f) / f;
    %% Calculate Electric Current Distribution in Cylindrical Coordinates
    J = zeros( size(E) );
    J(:, :, 1) = E(:, :, 2) .* ( cos( THe / 2 ) .^2 ) * ...
                 exp(-1j * k * 2 * f) / (Z * f);
    J(:, :, 2) = E(:, :, 3) .* ( cos( THe / 2 ) .^2 ) * ...
                 exp(-1j * k * 2 * f) / (Z * f);
end