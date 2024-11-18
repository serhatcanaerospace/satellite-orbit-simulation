function computePerturbations(a, e, i, RAAN, omega, mu)
    J2 = 1.08263e-3; % J2 katsay�s�
    Re = 6371; % D�nya'n�n yar��ap�

    % RAAN ve Arg�man De�i�im H�zlar�
    dRAAN_dt = -(3/2) * sqrt(mu) * J2 * (Re^2) * cos(i) / ((1 - e^2)^2 * a^(7/2));
    dOmega_dt = (3/4) * sqrt(mu) * J2 * (Re^2) * (5 * sin(i)^2 - 1) / ((1 - e^2)^2 * a^(7/2));

    % Zamanla De�i�imleri Hesaplama
    simulation_time = 24 * 3600; % 1 g�n (saniye)
    time_steps = linspace(0, simulation_time, 1000); % Zaman ad�mlar�
    RAAN_changes = RAAN + dRAAN_dt * time_steps;
    Omega_changes = omega + dOmega_dt * time_steps;

    % De�i�imleri �izdirme
    figure;
    plot(time_steps / 3600, rad2deg(RAAN_changes), 'r', 'LineWidth', 1.5); hold on;
    plot(time_steps / 3600, rad2deg(Omega_changes), 'b', 'LineWidth', 1.5);
    xlabel('Zaman (Saat)');
    ylabel('A��lar (Derece)');
    title('J2 Etkisi Alt�nda Y�r�nge De�i�imleri');
    legend('RAAN', 'Perigee Arg�man�');
    grid on;
end
