function computePerturbations(a, e, i, RAAN, omega, mu)
    J2 = 1.08263e-3; % J2 katsayýsý
    Re = 6371; % Dünya'nýn yarýçapý

    % RAAN ve Argüman Deðiþim Hýzlarý
    dRAAN_dt = -(3/2) * sqrt(mu) * J2 * (Re^2) * cos(i) / ((1 - e^2)^2 * a^(7/2));
    dOmega_dt = (3/4) * sqrt(mu) * J2 * (Re^2) * (5 * sin(i)^2 - 1) / ((1 - e^2)^2 * a^(7/2));

    % Zamanla Deðiþimleri Hesaplama
    simulation_time = 24 * 3600; % 1 gün (saniye)
    time_steps = linspace(0, simulation_time, 1000); % Zaman adýmlarý
    RAAN_changes = RAAN + dRAAN_dt * time_steps;
    Omega_changes = omega + dOmega_dt * time_steps;

    % Deðiþimleri Çizdirme
    figure;
    plot(time_steps / 3600, rad2deg(RAAN_changes), 'r', 'LineWidth', 1.5); hold on;
    plot(time_steps / 3600, rad2deg(Omega_changes), 'b', 'LineWidth', 1.5);
    xlabel('Zaman (Saat)');
    ylabel('Açýlar (Derece)');
    title('J2 Etkisi Altýnda Yörünge Deðiþimleri');
    legend('RAAN', 'Perigee Argümaný');
    grid on;
end
