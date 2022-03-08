function [ Jft ] = calculateCylFTCurrent( J, KX, KY, RHO, PH )
%calculateFTCurrent This function calculates Fourier Transform (FT) of
%current distribution in cylindrical coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
%   Note: optimize the fourier transform current distribution loop
    %% Calculate Element Spacing
    drho = RHO(1, 2) - RHO(1, 1);
    dph = PH(2, 1) - PH(1, 1);
    %% Calculate Size of Matrix Partition
    N = 1;
    m = 100;
    if size(KX, 1) > m
        N = size(KX, 1) / m;
        remN = rem( size(KX, 1), m );
        while remN ~= 0
            m = m + 1;
            N = size(KX, 1) / m;
            remN = rem( size(KX, 1), m );
        end
    end
    %% Alocate Space for Current Distribution Fourier Transform (FT)
    Jft = zeros([ size(KX), 3 ]);
    %% Transform Matricies
    KX = permute(KX, [3 4 1 2]);
    KY = permute(KY, [3 4 1 2]);
    %% Calculate Current Distribution Fourier Transform (FT)
    for i = 1 : N
        for n = 1 : N
            Jrho = sum( sum( J(:, :, 1) .* exp( 1j * ...
                 KX(:, :, (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m) .* ...
                 RHO .* cos(PH) ) .* exp( 1j * ...
                 KY(:, :, (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m) .* ...
                 RHO .* sin(PH) ) .* RHO ) ) * drho * dph;
            Jft( (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m, 1 ) = ...
                                                  permute(Jrho, [3 4 1 2]);
            Jph = sum( sum( J(:, :, 2) .* exp( 1j * ...
                 KX(:, :, (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m) .* ...
                 RHO .* cos(PH) ) .* exp( 1j * ...
                 KY(:, :, (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m) .* ...
                 RHO .* sin(PH) ) .* RHO ) ) * drho * dph;
            Jft( (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m, 2 ) = ...
                                                   permute(Jph, [3 4 1 2]);
            Jz = sum( sum( J(:, :, 3) .* exp( 1j * ...
                 KX(:, :, (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m) .* ...
                 RHO .* cos(PH) ) .* exp( 1j * ...
                 KY(:, :, (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m) .* ...
                 RHO .* sin(PH) ) .* RHO ) ) * drho * dph;
            Jft( (i - 1) * m + 1 : i * m, (n - 1) * m + 1 : n * m, 3 ) = ...
                                                    permute(Jz, [3 4 1 2]);
        end
    end
end