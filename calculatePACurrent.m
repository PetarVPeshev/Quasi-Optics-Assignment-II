function [ J, M ] = calculatePACurrent( E, Z, k, R, TH, f )
%calculatePACurrent This function calculates the equivalent aperture
%current of parabolic antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Electric Field Components with no Spherical Wave Part
    E = E * 2 * pi * R / exp(-1j * k * R);
    %% Calculate Magnetic Current Distribution in Cylindrical Coordinates
    M = zeros( size(E) );
    M(:, :, 1) = ( - E(:, :, 2) ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / f;
    M(:, :, 2) = ( E(:, :, 3)  ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / f;
    %% Calculate Electric Current Distribution in Cylindrical Coordinates
    J = zeros( size(E) );
    J(:, :, 1) = ( E(:, :, 3) ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / (Z * f);
    J(:, :, 2) = ( E(:, :, 2) ) .* ...
                 ( cos( TH / 2 ) .^2 ) * exp(-1j * k * 2 * f) / (Z * f);
end