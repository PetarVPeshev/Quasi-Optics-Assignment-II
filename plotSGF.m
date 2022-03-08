function [ fig ] = plotSGF( SGF, k, kx, componentSGF )
%plotSGF This function plots the components of the Spectral's Green
%Function (SGF)
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Plot Spectral's Green Function (SGF)
    fig = figure();
    subplot(3, 3, 1);
    plot( kx(1, :) / k, real( SGF(1, :, 1, 1) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 1, 1) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{xx}^{' componentSGF '}']);
    subplot(3, 3, 2);
    plot( kx(1, :) / k, real( SGF(1, :, 1, 2) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 1, 2) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{xy}^{' componentSGF '}']);
    subplot(3, 3, 3);
    plot( kx(1, :) / k, real( SGF(1, :, 1, 3) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 1, 3) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{xz}^{' componentSGF '}']);
    subplot(3, 3, 4);
    plot( kx(1, :) / k, real( SGF(1, :, 2, 1) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 2, 1) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{yx}^{' componentSGF '}']);
    subplot(3, 3, 5);
    plot( kx(1, :) / k, real( SGF(1, :, 2, 2) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 2, 2) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{yy}^{' componentSGF '}']);
    subplot(3, 3, 6);
    plot( kx(1, :) / k, real( SGF(1, :, 2, 3) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 2, 3) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{yz}^{' componentSGF '}']);
    subplot(3, 3, 7);
    plot( kx(1, :) / k, real( SGF(1, :, 3, 1) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 3, 1) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{zx}^{' componentSGF '}']);
    subplot(3, 3, 8);
    plot( kx(1, :) / k, real( SGF(1, :, 3, 2) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 3, 2) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{zy}^{' componentSGF '}']);
    subplot(3, 3, 9);
    plot( kx(1, :) / k, real( SGF(1, :, 3, 3) ), 'LineWidth', 3.0 );
    hold on;
    plot( kx(1, :) / k, imag( SGF(1, :, 3, 3) ), 'LineWidth', 3.0 );
    grid on;
    xlabel('k_{x} / k');
    ylabel(['G_{zz}^{' componentSGF '}']);
end