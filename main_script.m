% main_script.m

% Sabitler
mu = 398600; % Yerçekimi parametresi (km^3/s^2)

% TLE Dosyasýný Okuma
filename = 'C:\Users\TOSHIBA\Desktop\UNIFIED ENGINEERING_files\last-30-days.txt'; % TLE dosyasýnýn tam yolu
tle_lines = readTLE(filename); % TLE satýrlarýný oku

% TLE'deki Ýlk Uydunun Satýrlarýný Ayýrma
line0 = tle_lines{1}; % Uydu adý
line1 = tle_lines{2}; % 1. TLE satýrý
line2 = tle_lines{3}; % 2. TLE satýrý

% Kontrol: TLE Verilerini Doðru Okuduk mu?
disp('Uydu Adý:');
disp(line0);
disp('1. Satýr:');
disp(line1);
disp('2. Satýr:');
disp(line2);

% Ýlk TLE'yi Ýþleme
[a, e, i, RAAN, omega, M] = processTLE(line1, line2, mu);

% Yörüngeyi Görselleþtirme
visualizeOrbit(a, e, i, RAAN, omega, M, mu);

% Yörünge Pertürbasyonlarýný Hesaplama
computePerturbations(a, e, i, RAAN, omega, mu);

% Uydu Takýmý Simülasyonu
simulateSatelliteSwarm(tle_lines, mu);
