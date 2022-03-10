function [ J, M ] = calculatePACurrent( E, Z, k, R, TH, PH, f )
%calculatePACurrent This function calculates the equivalent aperture
%current of parabolic antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Electric Field Components with no Spherical Wave Part
    E = E * 2 * pi * R / exp(-1j * k * R);
    %% Create z-Component Normal in Cylindrical Coordinates
%     z = zeros( size(E) );
%     z(:, :, 3) = 1;
    %% Calculate Magnetic Current Distribution in Cylindrical Coordinates
    M = zeros( size(E) );
%     M(:, :, 1) = - E(:, :, 3) .* ( cos( TH / 2 ) .^2 ) * ...
%                  exp(-1j * k * 2 * f) / f;
%     M(:, :, 2) = - E(:, :, 2) .* ( cos( TH / 2 ) .^2 ) * ...
%                  exp(-1j * k * 2 * f) / f;
%     M = cross( M, z, 3 );
%     M = convertCylToCar(M, PH);
    M(:, :, 1) = ( - E(:, :, 2) .* cos(PH) - E(:, :, 3) .* sin(PH) ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / f;
    M(:, :, 2) = ( - E(:, :, 2) .* sin(PH) + E(:, :, 3) .* cos(PH)  ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / f;
    %% Calculate Electric Current Distribution in Cylindrical Coordinates
%     J = cross( M, z, 3 ) / Z;
%     J = convertCylToCar(J, PH);
    J = zeros( size(E) );
    J(:, :, 1) = ( E(:, :, 3) .* cos(PH) - E(:, :, 2) .* sin(PH) ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / (Z * f);
    J(:, :, 2) = ( E(:, :, 3) .* sin(PH) + E(:, :, 2) .* cos(PH) ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / (Z * f);
end