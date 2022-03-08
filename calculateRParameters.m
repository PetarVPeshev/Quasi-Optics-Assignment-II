function [ Dm, D, G ] = calculateRParameters( D, wl, Te, Ae )
%calculateRParameters This function calculates the maximum directivity,
%directivity, and gain of a reflector antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Area of Antenna
    A = pi * (D / 2)^2;
    %% Calculate Maximum Possible Directivity
    Dm = 4 * pi * A / ( wl ^ 2 );
    %% Calculate Directivity
    D = Dm * Te;
    %% Calculate Gain
    G = Dm * Ae;
end