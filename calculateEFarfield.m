function [ E ] = calculateEFarfield( ej_SGF, J, k, R, TH, kz )
%calculateFarfield This function calculates the electric far-field of
%certain current distribution and Spectral Green's Function
%   Detailed explanation goes here
%   Note: fix in-line documentation
%   Note: use matrix multiplication instead of for loop
    %% Calculate Electric Far-Field
    E = zeros([ size(TH), 3 ]);
    for i = 1:3
        for n = 1:3
            E(:, :, i) = E(:, :, i) + 1j .* kz .* ej_SGF(:, :, i, n) .* ...
                         J(:, :, n) * exp(-1j * k * R) / (2 * pi * R);
        end
    end
end