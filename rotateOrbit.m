function [x_rot, y_rot, z_rot] = rotateOrbit(x, y, z, inclination, RAAN)
    % Eðim ve RAAN'a göre yörünge noktalarýný döndürme
    R_inc = [1, 0, 0;
             0, cos(inclination), -sin(inclination);
             0, sin(inclination), cos(inclination)];

    R_RAAN = [cos(RAAN), -sin(RAAN), 0;
              sin(RAAN), cos(RAAN), 0;
              0, 0, 1];

    R_total = R_RAAN * R_inc;

    coords = R_total * [x; y; z];

    x_rot = coords(1, :);
    y_rot = coords(2, :);
    z_rot = coords(3, :);
end
