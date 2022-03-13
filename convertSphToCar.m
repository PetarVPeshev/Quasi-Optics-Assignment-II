function [ Acar ] = convertSphToCar( A, TH, PH )
%convertSphToCar This function converts vector from spherical to cartesian
%coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Create Transformation Matrix
    TM = zeros([ size(TH), 3, 3 ]);
    % Radial distance R
    TM(:, :, 1, 1) = sin(TH) .* cos(PH);
    TM(:, :, 1, 2) = cos(TH) .* cos(PH);
    TM(:, :, 1, 3) = - sin(PH);
    % Inclanation angle Theta
    TM(:, :, 2, 1) = sin(TH) .* sin(PH);
    TM(:, :, 2, 2) = cos(TH) .* sin(PH);
    TM(:, :, 2, 3) = cos(PH);
    % Azimuth angle Phi
    TM(:, :, 3, 1) = cos(TH);
    TM(:, :, 3, 2) = - sin(TH);
    %% Reshape Matricies
    TM = permute(TM, [3 4 1 2]);
    A = permute(A, [3 4 1 2]);
    %% Convert to Spherical Coordinates
    Acar = pagemtimes(TM, A);
    %% Reshape Matrix
    Acar = permute(Acar, [3 4 1 2]);
end