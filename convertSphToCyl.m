function [ Acyl ] = convertSphToCyl( A, TH )
%convertSphToCyl This function converts spherical to cylindrical
%coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Create Transformation Matrix
    TM = zeros([ size(TH), 3, 3 ]);
    % x-Coordinate
    TM(:, :, 1, 1) = sin(TH);
    TM(:, :, 1, 2) = 0;
    TM(:, :, 1, 3) = cos(TH);
    % y-Coordinate
    TM(:, :, 2, 1) = 0;
    TM(:, :, 2, 2) = 1;
    TM(:, :, 2, 3) = 0;
    % z-Coordinate
    TM(:, :, 3, 1) = cos(TH);
    TM(:, :, 3, 2) = 0;
    TM(:, :, 3, 3) = - sin(TH);
    %% Reshape Matricies
    TM = permute(TM, [3 4 1 2]);
    A = permute(A, [3 4 1 2]);
    %% Convert to Cartesian Coordinates
    Acyl = pagemtimes(TM, A);
    %% Reshape Matrix
    Acyl = permute(Acyl, [3 4 1 2]);
end