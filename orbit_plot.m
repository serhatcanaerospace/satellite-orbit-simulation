% Sabitler
mu = 398600; % Yer�ekimi parametresi, km^3/s^2 (D�nya i�in)

% TLE'den Orbital Elementleri ��z�mleme
inclination = str2double(line2(9:16)); % E�im (derece)
RAAN = str2double(line2(18:25)); % ��k�� d���m� boylam� (derece)
eccentricity = str2double(['0.' line2(27:33)]); % Eksantriklik
arg_perigee = str2double(line2(35:42)); % Perigee a��s� (derece)
mean_anomaly = str2double(line2(44:51)); % Ortalama anomali (derece)
mean_motion = str2double(line2(53:63)); % Ortalama hareket (devir/g�n)

% Ortalama Hareketten Yar� B�y�k Eksen Hesab�
n = mean_motion * 2 * pi / (24 * 3600); % Radyan/saniye
a = (mu / n^2)^(1/3); % Yar� b�y�k eksen (km)

% 3 Boyutlu Y�r�nge Noktalar�n� Hesaplama
t = linspace(0, 2*pi, 500); % Y�r�nge boyunca zaman noktalar�
theta = t; % Ger�ek Anomali (basitle�tirilmi�)

r = a * (1 - eccentricity^2) ./ (1 + eccentricity * cos(theta)); % Yar��ap
x = r .* cos(theta); % X ekseni
y = r .* sin(theta); % Y ekseni
z = zeros(size(x)); % Z ekseni (d�z y�r�nge)

% 3B G�rselle�tirme
figure;
plot3(x, y, z, 'b', 'LineWidth', 2); % Y�r�nge �izgisi
hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % D�nya merkezi
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Uydu Y�r�ngesi');
grid on;
axis equal;
% D�nya Modelini Ekleme
R_earth = 6371; % D�nya'n�n yar��ap� (km)
[X_earth, Y_earth, Z_earth] = sphere(50); % K�re koordinatlar�
X_earth = X_earth * R_earth; % X ekseni
Y_earth = Y_earth * R_earth; % Y ekseni
Z_earth = Z_earth * R_earth; % Z ekseni

% D�nya'y� �izdirme
hold on;
surf(X_earth, Y_earth, Z_earth, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% Uydu Animasyonu
satellite = plot3(0, 0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow'); % Uydu noktas�
num_steps = length(theta); % Zaman ad�mlar�

for i = 1:num_steps
    % Uydu pozisyonunu g�ncelle
    satellite.XData = r(i) * cos(theta(i));
    satellite.YData = r(i) * sin(theta(i));
    satellite.ZData = 0; % Z ekseni (d�z y�r�nge)
    
    % Grafi�i g�ncelle
    pause(0.01); % K���k bir bekleme (animasyon h�z�)
    drawnow;
end
% Perigee ve Apogee Hesab�
perigee = a * (1 - eccentricity); % En yak�n mesafe (km)
apogee = a * (1 + eccentricity); % En uzak mesafe (km)

% Y�r�nge Periyodu
T = 2 * pi * sqrt(a^3 / mu); % Saniye

% Sonu�lar� G�ster
fprintf('Perigee (En Yak�n Mesafe): %.2f km\n', perigee);
fprintf('Apogee (En Uzak Mesafe): %.2f km\n', apogee);
fprintf('Y�r�nge Periyodu: %.2f dakika\n', T / 60);
% Sabitler
J2 = 1.08263e-3; % J2 katsay�s�
Re = 6371; % D�nya'n�n yar��ap� (km)

% Pert�rbasyonlardan Dolay� De�i�imler
dRAAN_dt = -(3/2) * sqrt(mu) * J2 * (Re^2) * cosd(inclination) / ((1 - eccentricity^2)^2 * a^(7/2)); % RAAN de�i�imi (radyan/saniye)
dOmega_dt = (3/4) * sqrt(mu) * J2 * (Re^2) * (5 * sind(inclination)^2 - 1) / ((1 - eccentricity^2)^2 * a^(7/2)); % Perigee a��s� de�i�imi (radyan/saniye)

% Zaman �zerinde De�i�imleri Hesapla
simulation_time = 24 * 3600; % 1 g�n (saniye)
time_steps = linspace(0, simulation_time, 1000); % Zaman ad�mlar�
RAAN = RAAN + dRAAN_dt * time_steps; % ��k�� d���m� boylam�
Omega = arg_perigee + dOmega_dt * time_steps; % Perigee arg�man�

% De�i�imleri G�rselle�tirme
figure;
plot(time_steps / 3600, RAAN * 180 / pi, 'r', 'LineWidth', 1.5); hold on;
plot(time_steps / 3600, Omega * 180 / pi, 'b', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('A��lar (Derece)');
title('J2 Etkisi Alt�nda Y�r�nge De�i�imleri');
legend('RAAN', 'Perigee Arg�man�');
grid on;
% Y�r�nge Transferi: Hohmann Transferi
r1 = 7000; % Ba�lang�� y�r�nge yar��ap� (LEO), km
r2 = 42164; % Hedef y�r�nge yar��ap� (GEO), km

% 1. Delta-v (LEO -> Transfer Y�r�ngesi)
v1 = sqrt(mu / r1); % Ba�lang�� y�r�ngesindeki h�z
v_transfer1 = sqrt(mu * (2/r1 - 1/(r1 + r2))); % Transfer y�r�ngesindeki h�z
delta_v1 = abs(v_transfer1 - v1);

% 2. Delta-v (Transfer Y�r�ngesi -> GEO)
v2 = sqrt(mu / r2); % Hedef y�r�ngedeki h�z
v_transfer2 = sqrt(mu * (2/r2 - 1/(r1 + r2))); % Transfer y�r�ngesindeki h�z
delta_v2 = abs(v2 - v_transfer2);

% Toplam Delta-v
total_delta_v = delta_v1 + delta_v2;

% Sonu�lar� G�ster
fprintf('1. Delta-v (LEO -> Transfer Y�r�ngesi): %.2f km/s\n', delta_v1);
fprintf('2. Delta-v (Transfer Y�r�ngesi -> GEO): %.2f km/s\n', delta_v2);
fprintf('Toplam Delta-v: %.2f km/s\n', total_delta_v);
% �ki Uydunun Y�r�ngeleri Aras� �arp��ma Analizi
r_sat1 = a; % �lk uydunun y�r�nge yar��ap� (km)
r_sat2 = 8000; % �kinci uydunun y�r�nge yar��ap� (km)

% Y�r�nge Noktalar�n� Hesaplama
theta_sat2 = linspace(0, 2*pi, 500); % �kinci uydu y�r�ngesi
x_sat2 = r_sat2 * cos(theta_sat2);
y_sat2 = r_sat2 * sin(theta_sat2);

% Mesafe Hesaplama
distances = sqrt((x - x_sat2).^2 + (y - y_sat2).^2); % �ki y�r�nge aras�ndaki mesafeler

% �arp��ma Riskini Belirleme
min_distance = min(distances); % En k���k mesafe
collision_threshold = 10; % �arp��ma s�n�r� (km)

if min_distance < collision_threshold
    fprintf('Uyar�: �arp��ma riski! Minimum mesafe: %.2f km\n', min_distance);
else
    fprintf('G�vende: Minimum mesafe: %.2f km\n', min_distance);
end

% �kinci Y�r�ngeyi �izdirme
hold on;
plot3(x_sat2, y_sat2, zeros(size(x_sat2)), 'g--', 'LineWidth', 1.5);
legend('Uydu 1 Y�r�ngesi', 'D�nya', 'Uydu 2 Y�r�ngesi');

% Atmosferik S�r�klenme Analizi
Cd = 2.2; % S�r�klenme katsay�s� (tipik bir de�er)
A = 10; % Uydunun y�zey alan� (m^2)
m = 500; % Uydunun k�tlesi (kg)
rho = 1e-9; % Atmosfer yo�unlu�u (kg/m^3) - LEO i�in

% Ba�lang�� Y�r�ngesi
r = r1; % Ba�lang�� yar��ap� (km)
v = sqrt(mu / r); % H�z (km/s)
% G�ne� Radyasyonu Bas�nc� Analizi
S = 1361; % G�ne� radyasyon sabiti (W/m^2)
c = 3e8; % I��k h�z� (m/s)
P = S / c; % Bas�n� (N/m^2)
Cr = 1.5; % Yans�tma katsay�s� (tipik bir de�er)
A = 10; % Uydu y�zey alan� (m^2)
m = 500; % Uydu k�tlesi (kg)

% G�ne� Radyasyonu Bas�nc�ndan Kaynaklanan �vme
radiation_acceleration = Cr * P * A / m; % �vme (m/s^2)

% Zaman Ad�mlar�
dt = 60; % Zaman ad�m� (s)
simulation_time = 24 * 3600; % 1 g�n (saniye)
num_steps = simulation_time / dt;

% Y�r�nge Parametrelerindeki De�i�im
semi_major_axis = a; % Ba�lang�� yar� b�y�k ekseni (km)
semi_major_axis_history = zeros(1, num_steps);

for i = 1:num_steps
    % �vmenin Y�r�nge Elemanlar�na Etkisi
    delta_a = radiation_acceleration * dt; % Yar� b�y�k eksende de�i�im (km)
    semi_major_axis = semi_major_axis + delta_a / 1000; % Yeni yar� b�y�k eksen
    semi_major_axis_history(i) = semi_major_axis;
end

% Sonu�lar� G�rselle�tirme
time = (1:num_steps) * dt / 3600; % Zaman (saat)
figure;
plot(time, semi_major_axis_history - a, 'b', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Yar� B�y�k Eksendeki De�i�im (km)');
title('G�ne� Radyasyonu Bas�nc� Etkisi');
grid on;

% Zaman Ad�mlar�
dt = 60; % Zaman ad�m� (saniye)
simulation_time = 24 * 3600; % 1 g�n (saniye)
num_steps = simulation_time / dt;

% Y�r�nge Yar��ap�ndaki De�i�imi Hesaplama
radius_history = zeros(1, num_steps);
for i = 1:num_steps
    % Atmosferik S�r�klenmeden Kaynaklanan H�z De�i�imi
    drag_acceleration = 0.5 * rho * Cd * A / m * (v * 1000)^2; % m/s^2
    delta_v = drag_acceleration * dt; % H�z de�i�imi (m/s)
    v = v - delta_v / 1000; % Yeni h�z (km/s)
    
    % Y�r�nge Yar��ap�ndaki De�i�im
    r = mu / v^2; % Yar��ap (km)
    radius_history(i) = r;
end

% Sonu�lar� G�rselle�tirme
time = (1:num_steps) * dt / 3600; % Zaman (saat)
figure;
plot(time, radius_history - r1, 'r', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Y�r�nge Yar��ap�ndaki De�i�im (km)');
title('Atmosferik S�r�klenme Etkisi');
grid on;


% Uydu Tak�m� Sim�lasyonu
num_satellites = 6; % Uydu say�s�
radii = linspace(r1, r1 + 200, num_satellites); % Uydular�n farkl� y�r�nge yar��aplar�
theta_offset = linspace(0, 2*pi, num_satellites); % Uydular aras�ndaki faz fark�

% Zaman Ayarlar�
num_steps = 500;
theta = linspace(0, 2*pi, num_steps); % Y�r�nge boyunca a��lar

% Uydu Y�r�ngelerini Hesaplama
figure;
hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % D�nya merkezi
[X_earth, Y_earth, Z_earth] = sphere(50); % D�nya k�resi
surf(X_earth * Re, Y_earth * Re, Z_earth * Re, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

for i = 1:num_satellites
    r = radii(i);
    x = r * cos(theta + theta_offset(i));
    y = r * sin(theta + theta_offset(i));
    z = zeros(size(x)); % 2D y�r�nge
    
    plot3(x, y, z, 'LineWidth', 1.5); % Her uydunun y�r�ngesi
end

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Uydu Tak�m� Sim�lasyonu');
grid on;
axis equal;


% Ba�lang�� De�erleri
a = r1; % Ba�lang�� yar� b�y�k ekseni (km)
e = eccentricity; % Ba�lang�� eksantrikli�i
i = deg2rad(inclination); % Ba�lang�� e�imi (radyan)
RAAN = deg2rad(0); % Ba�lang�� RAAN
omega = deg2rad(arg_perigee); % Ba�lang�� perigee arg�man�
M = deg2rad(mean_anomaly); % Ba�lang�� ortalama anomali

% Zaman Parametreleri
dt = 60; % Zaman ad�m� (saniye)
simulation_time = 24 * 3600; % 1 g�n (saniye)
num_steps = simulation_time / dt;

% Pert�rbasyon Sabitleri
J2 = 1.08263e-3; % J2 katsay�s�
Cd = 2.2; % Atmosferik s�r�klenme katsay�s�
Cr = 1.5; % G�ne� radyasyonu yans�tma katsay�s�
rho = 1e-9; % Atmosfer yo�unlu�u (kg/m^3)
P = 1361 / c; % G�ne� radyasyonu bas�nc� (N/m^2)

% Y�r�nge De�i�imlerini Hesaplama
a_history = zeros(1, num_steps);
e_history = zeros(1, num_steps);
RAAN_history = zeros(1, num_steps);
omega_history = zeros(1, num_steps);

for step = 1:num_steps
    % J2 Pert�rbasyonu
    dRAAN_dt = -(3/2) * sqrt(mu) * J2 * (Re^2) * cos(i) / ((1 - e^2)^2 * a^(7/2));
    domega_dt = (3/4) * sqrt(mu) * J2 * (Re^2) * (5 * sin(i)^2 - 1) / ((1 - e^2)^2 * a^(7/2));
    
    % Atmosferik S�r�klenme
    drag_acceleration = 0.5 * rho * Cd * A / m * (sqrt(mu / a) * 1000)^2; % m/s^2
    delta_a_drag = -drag_acceleration * dt / (2 * sqrt(mu / a)); % Yar� b�y�k eksen de�i�imi
    
    % G�ne� Radyasyonu Bas�nc�
    delta_a_radiation = Cr * P * A / m * dt; % Yar� b�y�k eksen de�i�imi (m)
    
    % Y�r�nge Elemanlar�n� G�ncelle
    RAAN = RAAN + dRAAN_dt * dt; % RAAN de�i�imi
    omega = omega + domega_dt * dt; % Perigee arg�man� de�i�imi
    a = a + (delta_a_drag + delta_a_radiation / 1000); % Yar� b�y�k eksen
    e = e * exp(-drag_acceleration * dt); % Eksantriklik azalmas�
    
    % Sonu�lar� Kaydet
    a_history(step) = a;
    e_history(step) = e;
    RAAN_history(step) = rad2deg(RAAN);
    omega_history(step) = rad2deg(omega);
end

% Sonu�lar� G�rselle�tirme
time = (1:num_steps) * dt / 3600; % Zaman (saat)

figure;
subplot(3, 1, 1);
plot(time, a_history - r1, 'b', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Yar� B�y�k Eksendeki De�i�im (km)');
title('Yar� B�y�k Eksende De�i�im');

subplot(3, 1, 2);
plot(time, e_history, 'r', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('Eksantriklik');
title('Eksantriklikte De�i�im');

subplot(3, 1, 3);
plot(time, RAAN_history, 'g', 'LineWidth', 1.5);
xlabel('Zaman (Saat)');
ylabel('RAAN (Derece)');
title('RAAN De�i�imi');
% �oklu TLE ��leme
filename = 'last-30-days.txt'; % TLE dosya ad�
fileID = fopen(filename, 'r');
data = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

tle_lines = data{1};
num_tle = floor(length(tle_lines) / 3); % Her uydu i�in 3 sat�r

% Uydular�n Y�r�ngelerini Hesaplama
figure;
hold on;
[X_earth, Y_earth, Z_earth] = sphere(50);
surf(X_earth * Re, Y_earth * Re, Z_earth * Re, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

for i = 1:num_tle
    % TLE Verilerini ��z�mleme
    line1 = tle_lines{3*i-1};
    line2 = tle_lines{3*i};
    inclination = str2double(line2(9:16)); % E�im (derece)
    RAAN = str2double(line2(18:25)); % ��k�� d���m� boylam� (derece)
    eccentricity = str2double(['0.' line2(27:33)]); % Eksantriklik
    arg_perigee = str2double(line2(35:42)); % Perigee a��s� (derece)
    mean_anomaly = str2double(line2(44:51)); % Ortalama anomali (derece)
    mean_motion = str2double(line2(53:63)); % Ortalama hareket (devir/g�n)
    
    % Ortalama Hareketten Yar� B�y�k Eksen
    n = mean_motion * 2 * pi / (24 * 3600); % Radyan/saniye
    a = (mu / n^2)^(1/3); % Yar� b�y�k eksen (km)
    
    % Y�r�nge �izimi
    t = linspace(0, 2*pi, 500);
    r = a * (1 - eccentricity^2) ./ (1 + eccentricity * cos(t));
    x = r .* cos(t);
    y = r .* sin(t);
    z = zeros(size(x)); % D�z bir y�r�nge (ilk model)
    
    plot3(x, y, z, 'LineWidth', 1.5);
end

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('TLE Verileri ile Uydu Y�r�ngeleri');
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

