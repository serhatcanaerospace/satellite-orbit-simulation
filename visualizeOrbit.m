function visualizeOrbit(a, e, i, RAAN, omega, M, mu)
    % Yörünge Noktalarýný Hesaplama
    t = linspace(0, 2*pi, 500); % Gerçek Anomali
    r = a * (1 - e^2) ./ (1 + e * cos(t));
    x_orb = r .* cos(t);
    y_orb = r .* sin(t);
    z_orb = zeros(size(x_orb)); % Düz bir yörünge

    % Dünya Modelini Çizdirme
    figure;
    R_earth = 6371; % Dünya yarýçapý
    [X_earth, Y_earth, Z_earth] = sphere(50);
    surf(X_earth * R_earth, Y_earth * R_earth, Z_earth * R_earth, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;

    % Uydu Yörüngesini Çizdirme
    plot3(x_orb, y_orb, z_orb, 'b', 'LineWidth', 2);
    plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red'); % Dünya merkezi
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Uydu Yörüngesi');
    grid on;
    axis equal;
end
