function [a, e, i, RAAN, omega, M] = processTLE(line1, line2, mu)
    % TLE Verilerinden Yörünge Elemanlarýný Çözümleme
    
    % Eðim (inclination)
    inclination = str2double(line2(9:16)); % Derece
    if isnan(inclination)
        error('Eðim (inclination) deðeri TLE verilerinden alýnamadý.');
    end
    
    % Çýkýþ Düðümü Boylamý (RAAN)
    RAAN = str2double(line2(18:25)); % Derece
    if isnan(RAAN)
        error('RAAN deðeri TLE verilerinden alýnamadý.');
    end
    
    % Eksantriklik (eccentricity)
    eccentricity = str2double(['0.' line2(27:33)]); % Ondalýk formata çevrilir
    if isnan(eccentricity)
        error('Eksantriklik (eccentricity) deðeri TLE verilerinden alýnamadý.');
    end
    
    % Perigee Argümaný (arg_perigee)
    arg_perigee = str2double(line2(35:42)); % Derece
    if isnan(arg_perigee)
        error('Perigee argümaný (arg_perigee) deðeri TLE verilerinden alýnamadý.');
    end
    
    % Ortalama Anomali (mean_anomaly)
    mean_anomaly = str2double(line2(44:51)); % Derece
    if isnan(mean_anomaly)
        error('Ortalama anomali (mean_anomaly) deðeri TLE verilerinden alýnamadý.');
    end
    
    % Ortalama Hareket (mean_motion)
    mean_motion = str2double(line2(53:63)); % Devir/gün
    if isnan(mean_motion)
        error('Ortalama hareket (mean_motion) deðeri TLE verilerinden alýnamadý.');
    end

    % Ortalama Hareketten Yarý Büyük Eksen (semi-major axis)
    n = mean_motion * 2 * pi / (24 * 3600); % Radyan/saniye
    a = (mu / n^2)^(1/3); % Yarý büyük eksen (km)

    % Çýkýþ Deðerlerini Atama
    i = deg2rad(inclination); % Radyan
    RAAN = deg2rad(RAAN); % Radyan
    omega = deg2rad(arg_perigee); % Radyan
    M = deg2rad(mean_anomaly); % Radyan
    e = eccentricity; % Eksantriklik deðeri zaten 0 ile 1 arasýnda
end
