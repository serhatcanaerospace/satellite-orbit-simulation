% Sabitler
mu = 398600; % Yerçekimi parametresi, km^3/s^2 (Dünya için)

% TLE'den Orbital Elementleri Çözümleme
inclination = str2double(line2(9:16)); % Eðim (derece)
RAAN = str2double(line2(18:25)); % Çýkýþ düðümü boylamý (derece)
eccentricity = str2double(['0.' line2(27:33)]); % Eksantriklik
arg_perigee = str2double(line2(35:42)); % Perigee açýsý (derece)
mean_anomaly = str2double(line2(44:51)); % Ortalama anomali (derece)
mean_motion = str2double(line2(53:63)); % Ortalama hareket (devir/gün)

% Ortalama Hareketten Yarý Büyük Eksen Hesabý
n = mean_motion * 2 * pi / (24 * 3600); % Radyan/saniye
a = (mu / n^2)^(1/3); % Yarý büyük eksen (km)

% 3 Boyutlu Yörünge Noktalarýný Hesaplama
t = linspace(0, 2*pi, 500); % Yörünge boyunca zaman noktalarý
theta = t; % Gerçek Anomali (basitleþtirilmiþ)

r = a * (1 - eccentricity^2) ./ (1 + eccentricity * cos(theta)); % Yarýçap
x = r .* cos(theta); % X ekseni
y = r .* sin(theta); % Y ekseni
z = zeros(size(x)); % Z ekseni (düz yörünge)

% 3B Görselleþtirme
figure;
plot3(x, y, z, 'b', 'LineWidth', 2); % Yörünge çizgisi
hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Dünya merkezi
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Uydu Yörüngesi');
grid on;
axis equal;
% Dünya Modelini Ekleme
R_earth = 6371; % Dünya'nýn yarýçapý (km)
[X_earth, Y_earth, Z_earth] = sphere(50); % Küre koordinatlarý
X_earth = X_earth * R_earth; % X ekseni
Y_earth = Y_earth * R_earth; % Y ekseni
Z_earth = Z_earth * R_earth; % Z ekseni

% Dünya'yý Çizdirme
hold on;
surf(X_earth, Y_earth, Z_earth, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% Uydu Animasyonu
satellite = plot3(0, 0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow'); % Uydu noktasý
num_steps = length(theta); % Zaman adýmlarý

for i = 1:num_steps
    % Uydu pozisyonunu güncelle
    satellite.XData = r(i) * cos(theta(i));
    satellite.YData = r(i) * sin(theta(i));
    satellite.ZData = 0; % Z ekseni (düz yörünge)
    
    % Grafiði güncelle
    pause(0.01); % Küçük bir bekleme (animasyon hýzý)
    drawnow;
end
% Perigee ve Apogee Hesabý
perigee = a * (1 - eccentricity); % En yakýn mesafe (km)
apogee = a * (1 + eccentricity); % En uzak mesafe (km)

% Yörünge Periyodu
T = 2 * pi * sqrt(a^3 / mu); % Saniye

% Sonuçlarý Göster
fprintf('Perigee (En Yakýn Mesafe): %.2f km\n', perigee);
fprintf('Apogee (En Uzak Mesafe): %.2f km\n', apogee);
fprintf('Yörünge Periyodu: %.2f dakika\n', T / 60);
% Sabitler
J2 = 1.08263e-3; % J2 katsayýsý
Re = 6371; % Dünya'nýn yarýçapý (km)

% Pertürbasyonlardan Dolayý Deðiþimler
dRAAN_dt = -(3/2) * sqrt(mu) * J2 * (Re^2) * cosd(inclination) / ((1 - eccentricity^2)^2 * a^(7/2)); % RAAN deðiþimi (radyan/saniye)
dOmega_dt = (3/4) * sqrt(mu) * J2 * (Re^2) * (5 * sind(inclination)^2 - 1) / ((1 - eccentricity^2)^2 * a^(7/2)); % Perigee açýsý deðiþimi (radyan/saniye)

% Zaman Üzerinde Deðiþimleri Hesapla
simulation_time = 24 * 3600; % 1 gün (saniye)
time_steps = linspace(0, simulation_time, 1000); % Zaman adýmlarý
RAAN = RAAN + dRAAN_dt * time_steps; % Çýkýþ düðümü boylamý
Omega = arg_perigee + dOmega_dt * time_steps; % Perigee argümaný

% Deðiþimleri Görselleþtirme
figure;
plot(time_steps / 3600, RAAN * 180 / pi, 'r', 'LineWidth', 1.5); hold on;
plot(time_steps / 3600, Omega * 180 / pi, 'b', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Açýlar (Derece)');
title('J2 Etkisi Altýnda Yörünge Deðiþimleri');
legend('RAAN', 'Perigee Argümaný');
grid on;
% Yörünge Transferi: Hohmann Transferi
r1 = 7000; % Baþlangýç yörünge yarýçapý (LEO), km
r2 = 42164; % Hedef yörünge yarýçapý (GEO), km

% 1. Delta-v (LEO -> Transfer Yörüngesi)
v1 = sqrt(mu / r1); % Baþlangýç yörüngesindeki hýz
v_transfer1 = sqrt(mu * (2/r1 - 1/(r1 + r2))); % Transfer yörüngesindeki hýz
delta_v1 = abs(v_transfer1 - v1);

% 2. Delta-v (Transfer Yörüngesi -> GEO)
v2 = sqrt(mu / r2); % Hedef yörüngedeki hýz
v_transfer2 = sqrt(mu * (2/r2 - 1/(r1 + r2))); % Transfer yörüngesindeki hýz
delta_v2 = abs(v2 - v_transfer2);

% Toplam Delta-v
total_delta_v = delta_v1 + delta_v2;

% Sonuçlarý Göster
fprintf('1. Delta-v (LEO -> Transfer Yörüngesi): %.2f km/s\n', delta_v1);
fprintf('2. Delta-v (Transfer Yörüngesi -> GEO): %.2f km/s\n', delta_v2);
fprintf('Toplam Delta-v: %.2f km/s\n', total_delta_v);
% Ýki Uydunun Yörüngeleri Arasý Çarpýþma Analizi
r_sat1 = a; % Ýlk uydunun yörünge yarýçapý (km)
r_sat2 = 8000; % Ýkinci uydunun yörünge yarýçapý (km)

% Yörünge Noktalarýný Hesaplama
theta_sat2 = linspace(0, 2*pi, 500); % Ýkinci uydu yörüngesi
x_sat2 = r_sat2 * cos(theta_sat2);
y_sat2 = r_sat2 * sin(theta_sat2);

% Mesafe Hesaplama
distances = sqrt((x - x_sat2).^2 + (y - y_sat2).^2); % Ýki yörünge arasýndaki mesafeler

% Çarpýþma Riskini Belirleme
min_distance = min(distances); % En küçük mesafe
collision_threshold = 10; % Çarpýþma sýnýrý (km)

if min_distance < collision_threshold
    fprintf('Uyarý: Çarpýþma riski! Minimum mesafe: %.2f km\n', min_distance);
else
    fprintf('Güvende: Minimum mesafe: %.2f km\n', min_distance);
end

% Ýkinci Yörüngeyi Çizdirme
hold on;
plot3(x_sat2, y_sat2, zeros(size(x_sat2)), 'g--', 'LineWidth', 1.5);
legend('Uydu 1 Yörüngesi', 'Dünya', 'Uydu 2 Yörüngesi');

% Atmosferik Sürüklenme Analizi
Cd = 2.2; % Sürüklenme katsayýsý (tipik bir deðer)
A = 10; % Uydunun yüzey alaný (m^2)
m = 500; % Uydunun kütlesi (kg)
rho = 1e-9; % Atmosfer yoðunluðu (kg/m^3) - LEO için

% Baþlangýç Yörüngesi
r = r1; % Baþlangýç yarýçapý (km)
v = sqrt(mu / r); % Hýz (km/s)
% Güneþ Radyasyonu Basýncý Analizi
S = 1361; % Güneþ radyasyon sabiti (W/m^2)
c = 3e8; % Iþýk hýzý (m/s)
P = S / c; % Basýnç (N/m^2)
Cr = 1.5; % Yansýtma katsayýsý (tipik bir deðer)
A = 10; % Uydu yüzey alaný (m^2)
m = 500; % Uydu kütlesi (kg)

% Güneþ Radyasyonu Basýncýndan Kaynaklanan Ývme
radiation_acceleration = Cr * P * A / m; % Ývme (m/s^2)

% Zaman Adýmlarý
dt = 60; % Zaman adýmý (s)
simulation_time = 24 * 3600; % 1 gün (saniye)
num_steps = simulation_time / dt;

% Yörünge Parametrelerindeki Deðiþim
semi_major_axis = a; % Baþlangýç yarý büyük ekseni (km)
semi_major_axis_history = zeros(1, num_steps);

for i = 1:num_steps
    % Ývmenin Yörünge Elemanlarýna Etkisi
    delta_a = radiation_acceleration * dt; % Yarý büyük eksende deðiþim (km)
    semi_major_axis = semi_major_axis + delta_a / 1000; % Yeni yarý büyük eksen
    semi_major_axis_history(i) = semi_major_axis;
end

% Sonuçlarý Görselleþtirme
time = (1:num_steps) * dt / 3600; % Zaman (saat)
figure;
plot(time, semi_major_axis_history - a, 'b', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Yarý Büyük Eksendeki Deðiþim (km)');
title('Güneþ Radyasyonu Basýncý Etkisi');
grid on;

% Zaman Adýmlarý
dt = 60; % Zaman adýmý (saniye)
simulation_time = 24 * 3600; % 1 gün (saniye)
num_steps = simulation_time / dt;

% Yörünge Yarýçapýndaki Deðiþimi Hesaplama
radius_history = zeros(1, num_steps);
for i = 1:num_steps
    % Atmosferik Sürüklenmeden Kaynaklanan Hýz Deðiþimi
    drag_acceleration = 0.5 * rho * Cd * A / m * (v * 1000)^2; % m/s^2
    delta_v = drag_acceleration * dt; % Hýz deðiþimi (m/s)
    v = v - delta_v / 1000; % Yeni hýz (km/s)
    
    % Yörünge Yarýçapýndaki Deðiþim
    r = mu / v^2; % Yarýçap (km)
    radius_history(i) = r;
end

% Sonuçlarý Görselleþtirme
time = (1:num_steps) * dt / 3600; % Zaman (saat)
figure;
plot(time, radius_history - r1, 'r', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Yörünge Yarýçapýndaki Deðiþim (km)');
title('Atmosferik Sürüklenme Etkisi');
grid on;


% Uydu Takýmý Simülasyonu
num_satellites = 6; % Uydu sayýsý
radii = linspace(r1, r1 + 200, num_satellites); % Uydularýn farklý yörünge yarýçaplarý
theta_offset = linspace(0, 2*pi, num_satellites); % Uydular arasýndaki faz farký

% Zaman Ayarlarý
num_steps = 500;
theta = linspace(0, 2*pi, num_steps); % Yörünge boyunca açýlar

% Uydu Yörüngelerini Hesaplama
figure;
hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Dünya merkezi
[X_earth, Y_earth, Z_earth] = sphere(50); % Dünya küresi
surf(X_earth * Re, Y_earth * Re, Z_earth * Re, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

for i = 1:num_satellites
    r = radii(i);
    x = r * cos(theta + theta_offset(i));
    y = r * sin(theta + theta_offset(i));
    z = zeros(size(x)); % 2D yörünge
    
    plot3(x, y, z, 'LineWidth', 1.5); % Her uydunun yörüngesi
end

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Uydu Takýmý Simülasyonu');
grid on;
axis equal;


% Baþlangýç Deðerleri
a = r1; % Baþlangýç yarý büyük ekseni (km)
e = eccentricity; % Baþlangýç eksantrikliði
i = deg2rad(inclination); % Baþlangýç eðimi (radyan)
RAAN = deg2rad(0); % Baþlangýç RAAN
omega = deg2rad(arg_perigee); % Baþlangýç perigee argümaný
M = deg2rad(mean_anomaly); % Baþlangýç ortalama anomali

% Zaman Parametreleri
dt = 60; % Zaman adýmý (saniye)
simulation_time = 24 * 3600; % 1 gün (saniye)
num_steps = simulation_time / dt;

% Pertürbasyon Sabitleri
J2 = 1.08263e-3; % J2 katsayýsý
Cd = 2.2; % Atmosferik sürüklenme katsayýsý
Cr = 1.5; % Güneþ radyasyonu yansýtma katsayýsý
rho = 1e-9; % Atmosfer yoðunluðu (kg/m^3)
P = 1361 / c; % Güneþ radyasyonu basýncý (N/m^2)

% Yörünge Deðiþimlerini Hesaplama
a_history = zeros(1, num_steps);
e_history = zeros(1, num_steps);
RAAN_history = zeros(1, num_steps);
omega_history = zeros(1, num_steps);

for step = 1:num_steps
    % J2 Pertürbasyonu
    dRAAN_dt = -(3/2) * sqrt(mu) * J2 * (Re^2) * cos(i) / ((1 - e^2)^2 * a^(7/2));
    domega_dt = (3/4) * sqrt(mu) * J2 * (Re^2) * (5 * sin(i)^2 - 1) / ((1 - e^2)^2 * a^(7/2));
    
    % Atmosferik Sürüklenme
    drag_acceleration = 0.5 * rho * Cd * A / m * (sqrt(mu / a) * 1000)^2; % m/s^2
    delta_a_drag = -drag_acceleration * dt / (2 * sqrt(mu / a)); % Yarý büyük eksen deðiþimi
    
    % Güneþ Radyasyonu Basýncý
    delta_a_radiation = Cr * P * A / m * dt; % Yarý büyük eksen deðiþimi (m)
    
    % Yörünge Elemanlarýný Güncelle
    RAAN = RAAN + dRAAN_dt * dt; % RAAN deðiþimi
    omega = omega + domega_dt * dt; % Perigee argümaný deðiþimi
    a = a + (delta_a_drag + delta_a_radiation / 1000); % Yarý büyük eksen
    e = e * exp(-drag_acceleration * dt); % Eksantriklik azalmasý
    
    % Sonuçlarý Kaydet
    a_history(step) = a;
    e_history(step) = e;
    RAAN_history(step) = rad2deg(RAAN);
    omega_history(step) = rad2deg(omega);
end

% Sonuçlarý Görselleþtirme
time = (1:num_steps) * dt / 3600; % Zaman (saat)

figure;
subplot(3, 1, 1);
plot(time, a_history - r1, 'b', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Yarý Büyük Eksendeki Deðiþim (km)');
title('Yarý Büyük Eksende Deðiþim');

subplot(3, 1, 2);
plot(time, e_history, 'r', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Eksantriklik');
title('Eksantriklikte Deðiþim');

subplot(3, 1, 3);
plot(time, RAAN_history, 'g', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('RAAN (Derece)');
title('RAAN Deðiþimi');
% Çoklu TLE Ýþleme
filename = 'last-30-days.txt'; % TLE dosya adý
fileID = fopen(filename, 'r');
data = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

tle_lines = data{1};
num_tle = floor(length(tle_lines) / 3); % Her uydu için 3 satýr

% Uydularýn Yörüngelerini Hesaplama
figure;
hold on;
[X_earth, Y_earth, Z_earth] = sphere(50);
surf(X_earth * Re, Y_earth * Re, Z_earth * Re, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

for i = 1:num_tle
    % TLE Verilerini Çözümleme
    line1 = tle_lines{3*i-1};
    line2 = tle_lines{3*i};
    inclination = str2double(line2(9:16)); % Eðim (derece)
    RAAN = str2double(line2(18:25)); % Çýkýþ düðümü boylamý (derece)
    eccentricity = str2double(['0.' line2(27:33)]); % Eksantriklik
    arg_perigee = str2double(line2(35:42)); % Perigee açýsý (derece)
    mean_anomaly = str2double(line2(44:51)); % Ortalama anomali (derece)
    mean_motion = str2double(line2(53:63)); % Ortalama hareket (devir/gün)
    
    % Ortalama Hareketten Yarý Büyük Eksen
    n = mean_motion * 2 * pi / (24 * 3600); % Radyan/saniye
    a = (mu / n^2)^(1/3); % Yarý büyük eksen (km)
    
    % Yörünge Çizimi
    t = linspace(0, 2*pi, 500);
    r = a * (1 - eccentricity^2) ./ (1 + eccentricity * cos(t));
    x = r .* cos(t);
    y = r .* sin(t);
    z = zeros(size(x)); % Düz bir yörünge (ilk model)
    
    plot3(x, y, z, 'LineWidth', 1.5);
end

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('TLE Verileri ile Uydu Yörüngeleri');
grid on;
axis equal;



function [a, e, i, RAAN, omega, M] = processTLE(line1, line2, mu)
    inclination = str2double(line2(9:16));
    RAAN = str2double(line2(18:25));
    eccentricity = str2double(['0.' line2(27:33)]);
    arg_perigee = str2double(line2(35:42));
    mean_anomaly = str2double(line2(44:51));
    mean_motion = str2double(line2(53:63));
    n = mean_motion * 2 * pi / (24 * 3600);
    a = (mu / n^2)^(1/3);
    i = deg2rad(inclination);
    omega = deg2rad(arg_perigee);
    M = deg2rad(mean_anomaly);
end

