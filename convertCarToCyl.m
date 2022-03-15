function [ Acyl ] = convertCarToCyl( A, PH )
%convertCarToCyl This function converts vector from cartesian to
%cylindrical coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Create Transformation Matrix
    TM = zeros([ size(PH), 3, 3 ]);
    % Radial distance R
    TM(:, :, 1, 1) = cos(PH);
    TM(:, :, 1, 2) = sin(PH);
    % Inclanation angle Theta
    TM(:, :, 2, 1) = - sin(PH);
    TM(:, :, 2, 2) = cos(PH);
    % Azimuth angle Phi
    TM(:, :, 3, 3) = 1;
    %% Reshape Matricies
    TM = permute(TM, [3 4 1 2]);
    A = permute(A, [3 4 1 2]);
    %% Convert to Spherical Coordinates
    Acyl = pagemtimes(TM, A);
    %% Reshape Matrix
    Acyl = permute(Acyl, [3 4 1 2]);
end