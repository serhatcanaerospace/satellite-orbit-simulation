% main_script.m

% Sabitler
mu = 398600; % Yer�ekimi parametresi (km^3/s^2)

% TLE Dosyas�n� Okuma
filename = 'C:\Users\TOSHIBA\Desktop\UNIFIED ENGINEERING_files\last-30-days.txt'; % TLE dosyas�n�n tam yolu
tle_lines = readTLE(filename); % TLE sat�rlar�n� oku

% TLE'deki �lk Uydunun Sat�rlar�n� Ay�rma
line0 = tle_lines{1}; % Uydu ad�
line1 = tle_lines{2}; % 1. TLE sat�r�
line2 = tle_lines{3}; % 2. TLE sat�r�

% Kontrol: TLE Verilerini Do�ru Okuduk mu?
disp('Uydu Ad�:');
disp(line0);
disp('1. Sat�r:');
disp(line1);
disp('2. Sat�r:');
disp(line2);

% �lk TLE'yi ��leme
[a, e, i, RAAN, omega, M] = processTLE(line1, line2, mu);

% Y�r�ngeyi G�rselle�tirme
visualizeOrbit(a, e, i, RAAN, omega, M, mu);

% Y�r�nge Pert�rbasyonlar�n� Hesaplama
computePerturbations(a, e, i, RAAN, omega, mu);

% Uydu Tak�m� Sim�lasyonu
simulateSatelliteSwarm(tle_lines, mu);
