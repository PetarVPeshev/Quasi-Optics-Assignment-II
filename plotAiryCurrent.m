function [ figTh ] = plotAiryCurrent( J, TH )
%plotAiryCurrent This function plots the current distribution of circular
%source
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Map Vector to -90 to 90 Degrees
    ind = find( TH(1, :) > (90 * pi / 180) );
    JM = sqrt( J(1, :, 1).^2 + J(1, :, 2).^2 + J(1, :, 3).^2 );
    JM = circshift( JM, min(ind) );
    th = TH(1, :) * 180 / pi;
    th( min(ind) - 1 : end ) = th( min(ind) - 1 : end ) - 180;
    th = circshift( th, min(ind) );
    %% Plot 1D With Respect to Theta
    figTh = figure();
    plot( th, 20 * log10( JM ) - ...
          max( 20 * log10( JM ) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('\theta [deg]');
    ylabel('|FT(J)| [dB]');
    xlim([-90 90]);
    ylim([-40 0]);
    xticks((-90 : 10 : 90));
end