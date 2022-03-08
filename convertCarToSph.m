function [ Asph ] = convertCarToSph( A, TH, PH )
%convertCarToSph This function converts vector from cartesian to spherical
%coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Create Transformation Matrix
    TM = zeros([ size(TH), 3, 3 ]);
    % Radial distance R
    TM(:, :, 1, 1) = sin(TH) .* cos(PH);
    TM(:, :, 1, 2) = sin(TH) .* sin(PH);
    TM(:, :, 1, 3) = cos(TH);
    % Inclanation angle Theta
    TM(:, :, 2, 1) = cos(TH) .* cos(PH);
    TM(:, :, 2, 2) = cos(TH) .* sin(PH);
    TM(:, :, 2, 3) = - sin(TH);
    % Azimuth angle Phi
    TM(:, :, 3, 1) = -sin(TH);
    TM(:, :, 3, 2) = cos(PH);
    %% Reshape Matricies
    TM = permute(TM, [3 4 1 2]);
    A = permute(A, [3 4 1 2]);
    %% Convert to Spherical Coordinates
    Asph = pagemtimes(TM, A);
    %% Reshape Matrix
    Asph = permute(Asph, [3 4 1 2]);
end