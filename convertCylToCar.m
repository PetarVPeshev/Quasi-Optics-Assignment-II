function [ Acar ] = convertCylToCar( A, PH )
%convertCylToCar This function converts vector from cylindrical to
%cartesian coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Create Transformation Matrix
    TM = zeros([ size(PH), 3, 3 ]);
    % x-Coordinate
    TM(:, :, 1, 1) = cos(PH);
    TM(:, :, 1, 2) = - sin(PH);
    TM(:, :, 1, 3) = 0;
    % y-Coordinate
    TM(:, :, 2, 1) = sin(PH);
    TM(:, :, 2, 2) = cos(PH);
    TM(:, :, 2, 3) = 0;
    % z-Coordinate
    TM(:, :, 3, 1) = 0;
    TM(:, :, 3, 2) = 0;
    TM(:, :, 3, 3) = 1;
    %% Reshape Matricies
    TM = permute(TM, [3 4 1 2]);
    A = permute(A, [3 4 1 2]);
    %% Convert to Cartesian Coordinates
    Acar = pagemtimes(TM, A);
    %% Reshape Matrix
    Acar = permute(Acar, [3 4 1 2]);
end