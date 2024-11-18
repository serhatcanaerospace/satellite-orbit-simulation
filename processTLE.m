function [a, e, i, RAAN, omega, M] = processTLE(line1, line2, mu)
    % TLE Verilerinden Y�r�nge Elemanlar�n� ��z�mleme
    
    % E�im (inclination)
    inclination = str2double(line2(9:16)); % Derece
    if isnan(inclination)
        error('E�im (inclination) de�eri TLE verilerinden al�namad�.');
    end
    
    % ��k�� D���m� Boylam� (RAAN)
    RAAN = str2double(line2(18:25)); % Derece
    if isnan(RAAN)
        error('RAAN de�eri TLE verilerinden al�namad�.');
    end
    
    % Eksantriklik (eccentricity)
    eccentricity = str2double(['0.' line2(27:33)]); % Ondal�k formata �evrilir
    if isnan(eccentricity)
        error('Eksantriklik (eccentricity) de�eri TLE verilerinden al�namad�.');
    end
    
    % Perigee Arg�man� (arg_perigee)
    arg_perigee = str2double(line2(35:42)); % Derece
    if isnan(arg_perigee)
        error('Perigee arg�man� (arg_perigee) de�eri TLE verilerinden al�namad�.');
    end
    
    % Ortalama Anomali (mean_anomaly)
    mean_anomaly = str2double(line2(44:51)); % Derece
    if isnan(mean_anomaly)
        error('Ortalama anomali (mean_anomaly) de�eri TLE verilerinden al�namad�.');
    end
    
    % Ortalama Hareket (mean_motion)
    mean_motion = str2double(line2(53:63)); % Devir/g�n
    if isnan(mean_motion)
        error('Ortalama hareket (mean_motion) de�eri TLE verilerinden al�namad�.');
    end

    % Ortalama Hareketten Yar� B�y�k Eksen (semi-major axis)
    n = mean_motion * 2 * pi / (24 * 3600); % Radyan/saniye
    a = (mu / n^2)^(1/3); % Yar� b�y�k eksen (km)

    % ��k�� De�erlerini Atama
    i = deg2rad(inclination); % Radyan
    RAAN = deg2rad(RAAN); % Radyan
    omega = deg2rad(arg_perigee); % Radyan
    M = deg2rad(mean_anomaly); % Radyan
    e = eccentricity; % Eksantriklik de�eri zaten 0 ile 1 aras�nda
end
