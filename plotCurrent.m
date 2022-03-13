function [ figSurX, figSurY ] = plotCurrent( J, RHO, PH, compText )
%plotCurrent This function plots the current distribution in cartesian
%coordinates
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Interpolate Data for Smoother Surface
    rho_inter = linspace( min( RHO(1, :) ), max( RHO(1, :) ), ...
                          size( RHO, 2 ) * 5);
    ph_inter = linspace( min( PH(:, 1) ), max( PH(:, 1) ), ...
                          size( PH, 1 ) * 5);
    [ RHO_inter, PH_inter ] = meshgrid(rho_inter, ph_inter);
    J_inter = zeros([ size( RHO_inter ), 2 ]);
    J_inter(:, :, 1) = interp2(RHO, PH, J(:, :, 1), ...
                               RHO_inter, PH_inter, 'spline');
    J_inter(:, :, 2) = interp2(RHO, PH, J(:, :, 2), ...
                               RHO_inter, PH_inter, 'spline');
    %% Convert to Cartesian Coordinates
    X = RHO_inter .* cos(PH_inter);
    Y = RHO_inter .* sin(PH_inter);
    %% Plot 2D in XY Coordinates
    figSurX = figure();
    surface( X, Y, 20 * log10( abs( J_inter(:, :, 1) ) ) - ...
             max( max( 20 * log10( abs( J_inter(:, :, 1) ) ) ) ), ...
             'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    caxis([-40, 0]);
    view(0, 90);
    xlabel('X');
    ylabel('Y');
    zlabel(['|' compText '_{x}| [dB]']);
    zlim([-40 0]);
    %% Plot 2D in XY Coordinates
    figSurY = figure();
    surface( X, Y, 20 * log10( abs( J_inter(:, :, 2) ) ) - ...
             max( max( 20 * log10( abs( J_inter(:, :, 2) ) ) ) ), ...
             'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    caxis([-40, 0]);
    view(0, 90);
    xlabel('X');
    ylabel('Y');
    zlabel(['|' compText '_{y}| [dB]']);
    zlim([-40 0]);
end