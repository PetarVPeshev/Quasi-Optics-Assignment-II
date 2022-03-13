function [ figSur ] = plotFarfield( F, TH, PH )
%plotFarfield This function plots the far-farfield in UV representation
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Extract Dimensions of Theta and Phi
    Nth = size(TH, 2);
    Nph = size(PH, 1);
    %% Interpolate Data for Smoother Surface
    th_inter = linspace(eps, pi, Nth * 5);
    ph_inter = linspace(eps, 2 * pi, Nph * 5);
    [ TH_inter, PH_inter ] = meshgrid(th_inter, ph_inter);
    Fr_inter = interp2(TH, PH, F(:, :, 1), TH_inter, PH_inter, 'spline');
    Fth_inter = interp2(TH, PH, F(:, :, 2), TH_inter, PH_inter, 'spline');
    Fphi_inter = interp2(TH, PH, F(:, :, 3), TH_inter, PH_inter, 'spline');
    F_inter = sqrt( abs( Fr_inter ).^2 + abs( Fth_inter ).^2 + ...
                    abs( Fphi_inter ).^2 );
    %% UV Representation Coordinates
    U = sin(TH_inter) .* cos(PH_inter);
    V = sin(TH_inter) .* sin(PH_inter);
    %% Plot UV Representation of Far-Field
    figSur = figure();
    surface(U, V, 20 * log10( abs( F_inter ) ) - ...
            max( max( 20 * log10( abs( F_inter ) ) ) ), 'LineStyle', 'none' );
    grid on;
    colormap('jet');
    colorbar;
    caxis([-10, 0]);
    view(0, 90);
    xlabel('U');
    ylabel('V');
    zlabel('|E| [dB]');
    zlim([-10 0]);
end