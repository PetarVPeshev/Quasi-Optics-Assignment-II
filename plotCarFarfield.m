function [ figSur ] = plotCarFarfield( F, RHO, PH )
%plotCarFarfield This function plots the far-field in Cylindrical coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Extract Dimension of RHO and Phi
    Nrho = size(RHO, 2);
    Nph = size(PH, 1);
    %% Interpolate Data for Smoother Surface
    rho_inter = linspace(eps, max( RHO, [], [1  2] ), Nrho * 5);
    ph_inter = linspace(eps, 2 * pi, Nph * 5);
    [ RHO_inter, PH_inter ] = meshgrid(rho_inter, ph_inter);
    Frho_inter = interp2(RHO, PH, F(:, :, 1), RHO_inter, PH_inter, 'spline');
    Fph_inter = interp2(RHO, PH, F(:, :, 2), RHO_inter, PH_inter, 'spline');
    Fz_inter = interp2(RHO, PH, F(:, :, 3), RHO_inter, PH_inter, 'spline');
    F_inter = sqrt( abs( Frho_inter ).^2 + abs( Fph_inter ).^2 + ...
                    abs( Fz_inter ).^2 );
    %% Calculate X and Y Coordinates
    X = RHO_inter .* cos(PH_inter);
    Y = RHO_inter .* sin(PH_inter);
    %% Plot Far-Field in XY Coordinates
    figSur = figure();
    surface(X, Y, 20 * log10( abs( F_inter ) ) - ...
            max( max( 20 * log10( abs( F_inter ) ) ) ), 'LineStyle', 'none' );
    grid on;
    colormap('jet');
    colorbar;
    caxis([-10, 0]);
    view(0, 90);
    xlabel('X');
    ylabel('Y');
    zlabel('|E| [dB]');
    zlim([-10 0]);
end