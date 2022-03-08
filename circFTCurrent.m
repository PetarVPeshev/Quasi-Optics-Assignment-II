function [ J, Ja ] = circFTCurrent( k, A, TH, r, p )
%circFTCurrent This function computes the current distribution vector for
%uniform circularly distributed field
%   Detailed explanation goes here
%   Note: Fix in-line documentation
    %% Default Polarization in x
    switch nargin
        case 4
            p = [1 0 0];
        otherwise
            if sqrt( sum( p.^2 ) ) ~= 1
                error('Error. Invalid orientation vector.');
            end
    end
    %% Calculate Current Distribution Magnitude
    Ja = A * 2 * pi * r^2 * besselj(1, k * r * sin(TH)) ./ ...
         ( k * r * sin(TH) );
    %% Calculate Current Distribution Vector
    J = zeros([ size( TH ), 3 ]);
    for i = 1:3
        J(:, :, i) = Ja * p(i);
    end
end