function simulateSatelliteSwarm(tle_lines, mu)
    num_tle = floor(length(tle_lines) / 3); % Uydu sayýsý
    figure;
    hold on;

    % Dünya'yý Çizdirme
    R_earth = 6371; % Dünya yarýçapý
    [X_earth, Y_earth, Z_earth] = sphere(50);
    surf(X_earth * R_earth, Y_earth * R_earth, Z_earth * R_earth, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % Her Uydu Ýçin Yörünge Çizimi
    for i = 1:num_tle
        line1 = tle_lines{3*i - 1};
        line2 = tle_lines{3*i};
        [a, e, ~, RAAN, omega, ~] = processTLE(line1, line2, mu);

        t = linspace(0, 2*pi, 200);
        r = a * (1 - e^2) ./ (1 + e * cos(t));
        x = r .* cos(t);
        y = r .* sin(t);
        z = zeros(size(x)); % Düz yörünge

        plot3(x, y, z, 'LineWidth', 1.5);
    end

    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Uydu Takýmý Simülasyonu');
    grid on;
    axis equal;
end
