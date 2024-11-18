function visualizeOrbit(a, e, i, RAAN, omega, M, mu)
    % Y�r�nge Noktalar�n� Hesaplama
    t = linspace(0, 2*pi, 500); % Ger�ek Anomali
    r = a * (1 - e^2) ./ (1 + e * cos(t));
    x_orb = r .* cos(t);
    y_orb = r .* sin(t);
    z_orb = zeros(size(x_orb)); % D�z bir y�r�nge

    % D�nya Modelini �izdirme
    figure;
    R_earth = 6371; % D�nya yar��ap�
    [X_earth, Y_earth, Z_earth] = sphere(50);
    surf(X_earth * R_earth, Y_earth * R_earth, Z_earth * R_earth, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;

    % Uydu Y�r�ngesini �izdirme
    plot3(x_orb, y_orb, z_orb, 'b', 'LineWidth', 2);
    plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red'); % D�nya merkezi
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Uydu Y�r�ngesi');
    grid on;
    axis equal;
end
