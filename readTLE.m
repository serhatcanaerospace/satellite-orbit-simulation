function tle_lines = readTLE(filename)
    % TLE dosyasýný okuma
    fileID = fopen(filename, 'r'); % Dosyayý aç
    if fileID == -1
        error('TLE dosyasý açýlamadý. Yolu kontrol edin: %s', filename);
    end
    data = textscan(fileID, '%s', 'Delimiter', '\n'); % Satýrlarý oku
    fclose(fileID); % Dosyayý kapat

    tle_lines = data{1}; % Satýrlarý hücre dizisine aktar
end
