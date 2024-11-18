function tle_lines = readTLE(filename)
    % TLE dosyas�n� okuma
    fileID = fopen(filename, 'r'); % Dosyay� a�
    if fileID == -1
        error('TLE dosyas� a��lamad�. Yolu kontrol edin: %s', filename);
    end
    data = textscan(fileID, '%s', 'Delimiter', '\n'); % Sat�rlar� oku
    fclose(fileID); % Dosyay� kapat

    tle_lines = data{1}; % Sat�rlar� h�cre dizisine aktar
end
