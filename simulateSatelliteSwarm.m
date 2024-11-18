function simulateSatelliteSwarm(tle_lines, mu)
    num_tle = floor(length(tle_lines) / 3); % Uydu say�s�
    figure;
    hold on;

    % D�nya'y� �izdirme
    R_earth = 6371; % D�nya yar��ap�
    [X_earth, Y_earth, Z_earth] = sphere(50);
    surf(X_earth * R_earth, Y_earth * R_earth, Z_earth * R_earth, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % Her Uydu ��in Y�r�nge �izimi
    for i = 1:num_tle
        line1 = tle_lines{3*i - 1};
        line2 = tle_lines{3*i};
        [a, e, ~, RAAN, omega, ~] = processTLE(line1, line2, mu);

        t = linspace(0, 2*pi, 200);
        r = a * (1 - e^2) ./ (1 + e * cos(t));
        x = r .* cos(t);
        y = r .* sin(t);
        z = zeros(size(x)); % D�z y�r�nge

        plot3(x, y, z, 'LineWidth', 1.5);
    end

    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Uydu Tak�m� Sim�lasyonu');
    grid on;
    axis equal;
end
