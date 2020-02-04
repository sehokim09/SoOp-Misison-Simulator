clear all;

Tx = 'Beidou'; % MUOS, ORBCOMM, GPS, Glonass, Galileo, Beidou, NOAA, METEOR
Src = 'Web'; % TLE source - Web / File

if(strcmp(Tx,'MUOS'))
    %% MUOS
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_MUOS.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://www.celestrak.com/NORAD/elements/geo.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_MUOS.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_MUOS.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');
    
    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);
        ID = sscanf(Label,'MUOS-%d');
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);
        
        if(isempty(ID))
        else
            % TLE input
            Label = sprintf('MUOS %d                  \n',ID);        
            fprintf(fileID_TLE,'%s',Label);
            fprintf(fileID_TLE,'%s',TLE1);
            fprintf(fileID_TLE,'%s',TLE2);

            % Orbit File
            Orb_Filename = sprintf('Orb_MUOS_%d.txt',ID);
            system(['cp Orb.txt ' Orb_Filename]);

            s1 = sprintf('MUOS %d	              !  Description\n', ID);
            s2 = sprintf('TLE  "MUOS %d"          !  TLE or TRV format, Label to find in file\n', ID);
            s3 = sprintf('"Inp_TLE_MUOS.txt"               !  File name\n');
            OverwriteLineInFile(Orb_Filename,2,s1);
            OverwriteLineInFile(Orb_Filename,23,s2);
            OverwriteLineInFile(Orb_Filename,24,s3);

            % SC File
            SC_Filename = sprintf('SC_MUOS_%d.txt',ID);
            system(['cp SC.txt ' SC_Filename]);

            s1 = sprintf('MUOS %d	              !  Description\n', ID);
            s2 = sprintf('"MUOS%d"                        !  Label\n', ID);
            s3 = sprintf('GenTx1SpriteAlpha.ppm         !  Sprite File Name\n');

            OverwriteLineInFile(SC_Filename,2,s1);
            OverwriteLineInFile(SC_Filename,3,s2);
            OverwriteLineInFile(SC_Filename,4,s3);
        end
    end
    
elseif(strcmp(Tx,'ORBCOMM'))
    %% ORBCOMM
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_ORBCOMM.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://celestrak.com/NORAD/elements/orbcomm.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_ORBCOMM.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_ORBCOMM.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');
    
    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);
        ID = sscanf(Label,'ORBCOMM FM%d');
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);
        
        if(isempty(ID))
        elseif(ismember(ID,[1,2,3,16,17,22,24,25,26,28,29,33,37,38,39,40,41,101,104,105,106,111,119])) % No longer operational
        else
            % TLE input
            Label = sprintf('ORBCOMM FM %d                 \n',ID);        
            fprintf(fileID_TLE,'%s',Label);
            fprintf(fileID_TLE,'%s',TLE1);
            fprintf(fileID_TLE,'%s',TLE2);

            % Orbit File
            Orb_Filename = sprintf('Orb_ORBCOMM_FM%d.txt',ID);
            system(['cp Orb.txt ' Orb_Filename]);

            s1 = sprintf('ORBCOMM FM %d	              !  Description\n', ID);
            s2 = sprintf('TLE  "ORBCOMM FM %d"          !  TLE or TRV format, Label to find in file\n', ID);
            s3 = sprintf('"Inp_TLE_ORBCOMM.txt"               !  File name\n');
            OverwriteLineInFile(Orb_Filename,2,s1);
            OverwriteLineInFile(Orb_Filename,23,s2);
            OverwriteLineInFile(Orb_Filename,24,s3);

            % SC File
            SC_Filename = sprintf('SC_ORBCOMM_FM%d.txt',ID);
            system(['cp SC.txt ' SC_Filename]);

            s1 = sprintf('ORBCOMM FM %d	              !  Description\n', ID);
            s2 = sprintf('"FM%d"                        !  Label\n', ID);
            s3 = sprintf('GenTx2SpriteAlpha.ppm         !  Sprite File Name\n');

            OverwriteLineInFile(SC_Filename,2,s1);
            OverwriteLineInFile(SC_Filename,3,s2);
            OverwriteLineInFile(SC_Filename,4,s3);
        end
    end
    
elseif(strcmp(Tx,'GPS'))
    %% GPS  
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_GPS.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://celestrak.com/NORAD/elements/gps-ops.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_GPS.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_GPS.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');

    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);
        k = strfind(Label,'PRN');
        ID = sscanf(Label(k:end),'PRN %d');
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);

        % TLE input
        Label = sprintf('GPS PRN %02d               \n',ID);        
        fprintf(fileID_TLE,'%s',Label);
        fprintf(fileID_TLE,'%s',TLE1);
        fprintf(fileID_TLE,'%s',TLE2);
        
        % Orbit File
        Orb_Filename = sprintf('Orb_GPS_PRN%d.txt',ID);
        system(['cp Orb.txt ' Orb_Filename]);

        s1 = sprintf('GPS PRN %d	              !  Description\n', ID);
        s2 = sprintf('TLE  "GPS PRN %02d"          !  TLE or TRV format, Label to find in file\n', ID);
        s3 = sprintf('"Inp_TLE_GPS.txt"               !  File name\n');
        OverwriteLineInFile(Orb_Filename,2,s1);
        OverwriteLineInFile(Orb_Filename,23,s2);
        OverwriteLineInFile(Orb_Filename,24,s3);

        % SC File
        SC_Filename = sprintf('SC_GPS_PRN%d.txt',ID);
        system(['cp SC.txt ' SC_Filename]);

        s1 = sprintf('GPS PRN %d	              !  Description\n', ID);
        s2 = sprintf('"GPS%d"                        !  Label\n', ID);
        s3 = sprintf('GenTx3SpriteAlpha.ppm         !  Sprite File Name\n');

        OverwriteLineInFile(SC_Filename,2,s1);
        OverwriteLineInFile(SC_Filename,3,s2);
        OverwriteLineInFile(SC_Filename,4,s3);
    end
    
elseif(strcmp(Tx,'Glonass')) 
    %% Glonass
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_Glonass.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://celestrak.com/NORAD/elements/glo-ops.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_Glonass.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_Glonass.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');

    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);        
        ID = sscanf(Label,'COSMOS %04d');
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);

        % TLE input
        Label = sprintf('Glonass COSMOS %04d             \n',ID);        
        fprintf(fileID_TLE,'%s',Label);
        fprintf(fileID_TLE,'%s',TLE1);
        fprintf(fileID_TLE,'%s',TLE2);
        
        % Orbit File
        Orb_Filename = sprintf('Orb_Glonass_COSMOS%04d.txt',ID);
        system(['cp Orb.txt ' Orb_Filename]);

        s1 = sprintf('Glonass COSMOS %04d	              !  Description\n', ID);
        s2 = sprintf('TLE  "Glonass COSMOS %04d"          !  TLE or TRV format, Label to find in file\n', ID);
        s3 = sprintf('"Inp_TLE_Glonass.txt"               !  File name\n');
        OverwriteLineInFile(Orb_Filename,2,s1);
        OverwriteLineInFile(Orb_Filename,23,s2);
        OverwriteLineInFile(Orb_Filename,24,s3);

        % SC File
        SC_Filename = sprintf('SC_Glonass_COSMOS%04d.txt',ID);
        system(['cp SC.txt ' SC_Filename]);

        s1 = sprintf('Glonass COSMOS %04d	              !  Description\n', ID);
        s2 = sprintf('"COSMOS%d"                        !  Label\n', ID);
        s3 = sprintf('GenTx4SpriteAlpha.ppm         !  Sprite File Name\n');

        OverwriteLineInFile(SC_Filename,2,s1);
        OverwriteLineInFile(SC_Filename,3,s2);
        OverwriteLineInFile(SC_Filename,4,s3);
    end

elseif(strcmp(Tx,'Galileo')) 
    %% Galileo    
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_Galileo.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://celestrak.com/NORAD/elements/galileo.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_Galileo.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_Galileo.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');

    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);
        k = strfind(Label,'PRN');
        ID = sscanf(Label(k:end),'PRN E%d)');
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);
        
        if(isempty(ID))
        else
            % TLE input
            Label = sprintf('Galileo E%02d             \n',ID);        
            fprintf(fileID_TLE,'%s',Label);
            fprintf(fileID_TLE,'%s',TLE1);
            fprintf(fileID_TLE,'%s',TLE2);

            % Orbit File
            Orb_Filename = sprintf('Orb_Galileo_E%d.txt',ID);
            system(['cp Orb.txt ' Orb_Filename]);

            s1 = sprintf('Galileo E%d	              !  Description\n', ID);
            s2 = sprintf('TLE  "Galileo E%02d"          !  TLE or TRV format, Label to find in file\n', ID);
            s3 = sprintf('"Inp_TLE_Galileo.txt"               !  File name\n');
            OverwriteLineInFile(Orb_Filename,2,s1);
            OverwriteLineInFile(Orb_Filename,23,s2);
            OverwriteLineInFile(Orb_Filename,24,s3);

            % SC File
            SC_Filename = sprintf('SC_Galileo_E%d.txt',ID);
            system(['cp SC.txt ' SC_Filename]);

            s1 = sprintf('Galileo E%d	              !  Description\n', ID);
            s2 = sprintf('"E%d"                        !  Label\n', ID);
            s3 = sprintf('GenTx5SpriteAlpha.ppm         !  Sprite File Name\n');

            OverwriteLineInFile(SC_Filename,2,s1);
            OverwriteLineInFile(SC_Filename,3,s2);
            OverwriteLineInFile(SC_Filename,4,s3);
        end
    end
    
elseif(strcmp(Tx,'Beidou')) 
    %% Beidou   
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_Beidou.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://celestrak.com/NORAD/elements/beidou.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_Beidou.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_Beidou.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');

    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);
        k = strfind(Label,'(C');
        ID = sscanf(Label(k:end),'(C%d)');
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);
        
        if(isempty(ID))
        else
            % TLE input
            Label = sprintf('Beidou C%02d              \n',ID);        
            fprintf(fileID_TLE,'%s',Label);
            fprintf(fileID_TLE,'%s',TLE1);
            fprintf(fileID_TLE,'%s',TLE2);

            % Orbit File
            Orb_Filename = sprintf('Orb_Beidou_C%d.txt',ID);
            system(['cp Orb.txt ' Orb_Filename]);

            s1 = sprintf('Beidou C%d	              !  Description\n', ID);
            s2 = sprintf('TLE  "Beidou C%02d"          !  TLE or TRV format, Label to find in file\n', ID);
            s3 = sprintf('"Inp_TLE_Beidou.txt"               !  File name\n');
            OverwriteLineInFile(Orb_Filename,2,s1);
            OverwriteLineInFile(Orb_Filename,23,s2);
            OverwriteLineInFile(Orb_Filename,24,s3);

            % SC File
            SC_Filename = sprintf('SC_Beidou_C%d.txt',ID);
            system(['cp SC.txt ' SC_Filename]);

            s1 = sprintf('Beidou C%d	              !  Description\n', ID);
            s2 = sprintf('"C%d"                        !  Label\n', ID);
            s3 = sprintf('GenTx6SpriteAlpha.ppm         !  Sprite File Name\n');

            OverwriteLineInFile(SC_Filename,2,s1);
            OverwriteLineInFile(SC_Filename,3,s2);
            OverwriteLineInFile(SC_Filename,4,s3);
        end
    end
elseif(strcmp(Tx,'NOAA')) 
    %% NOAA 
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_NOAA.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://celestrak.com/NORAD/elements/noaa.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_NOAA.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_NOAA.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');

    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);
        ID = sscanf(Label,'NOAA %d');
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);
        
        if(isempty(ID))
        else
            % TLE input
            Label = sprintf('NOAA %02d                  \n',ID);        
            fprintf(fileID_TLE,'%s',Label);
            fprintf(fileID_TLE,'%s',TLE1);
            fprintf(fileID_TLE,'%s',TLE2);

            % Orbit File
            Orb_Filename = sprintf('Orb_NOAA_%d.txt',ID);
            system(['cp Orb.txt ' Orb_Filename]);

            s1 = sprintf('NOAA %d	              !  Description\n', ID);
            s2 = sprintf('TLE  "NOAA %02d"          !  TLE or TRV format, Label to find in file\n', ID);
            s3 = sprintf('"Inp_TLE_NOAA.txt"               !  File name\n');
            OverwriteLineInFile(Orb_Filename,2,s1);
            OverwriteLineInFile(Orb_Filename,23,s2);
            OverwriteLineInFile(Orb_Filename,24,s3);

            % SC File
            SC_Filename = sprintf('SC_NOAA_%d.txt',ID);
            system(['cp SC.txt ' SC_Filename]);

            s1 = sprintf('NOAA %d	              !  Description\n', ID);
            s2 = sprintf('"NOAA%d"                        !  Label\n', ID);
            s3 = sprintf('GenTx7SpriteAlpha.ppm         !  Sprite File Name\n');

            OverwriteLineInFile(SC_Filename,2,s1);
            OverwriteLineInFile(SC_Filename,3,s2);
            OverwriteLineInFile(SC_Filename,4,s3);
        end
    end
elseif(strcmp(Tx,'METEOR')) 
    %% METEOR
    if(strcmp(Src,'Web'))
        src_Filename = 'TLE_from_CelesTrak_METEOR.txt';
        fileID_src = fopen(src_Filename, 'w');
        data = webread('https://celestrak.com/NORAD/elements/weather.txt');
        fprintf(fileID_src, '%s', data);
        fclose(fileID_src);
    end
    input_Filename = 'TLE_from_CelesTrak_METEOR.txt';
    fileID_input = fopen(input_Filename, 'r');
    
    TLE_Filename = 'Inp_TLE_METEOR.txt';
    fileID_TLE = fopen(TLE_Filename, 'w');

    while ~feof(fileID_input)
        Label = fgets(fileID_input,26);
        TLE1 = fgets(fileID_input,512);
        TLE2 = fgets(fileID_input,512);
        
        if(strcmp(TLE1(3:7),'35865'))
            flag = 1;
            ID = 'M1';
        elseif(strcmp(TLE1(3:7),'40069'))
            flag = 1;
            ID = 'M2';
        elseif(strcmp(TLE1(3:7),'44387'))
            flag = 1;
            ID = 'M2-2';
        else
            flag = 0;
        end
        
        if(flag)
            % TLE input
            Label = sprintf('METEOR-%s                \n',ID);        
            fprintf(fileID_TLE,'%s',Label);
            fprintf(fileID_TLE,'%s',TLE1);
            fprintf(fileID_TLE,'%s',TLE2);

            % Orbit File
            Orb_Filename = sprintf('Orb_METEOR_%s.txt',ID);
            system(['cp Orb.txt ' Orb_Filename]);

            s1 = sprintf('METEOR %s	              !  Description\n', ID);
            s2 = sprintf('TLE  "METEOR-%s"          !  TLE or TRV format, Label to find in file\n', ID);
            s3 = sprintf('"Inp_TLE_METEOR.txt"               !  File name\n');
            OverwriteLineInFile(Orb_Filename,2,s1);
            OverwriteLineInFile(Orb_Filename,23,s2);
            OverwriteLineInFile(Orb_Filename,24,s3);

            % SC File
            SC_Filename = sprintf('SC_METEOR_%s.txt',ID);
            system(['cp SC.txt ' SC_Filename]);

            s1 = sprintf('METEOR %s	              !  Description\n', ID);
            s2 = sprintf('"%s"                        !  Label\n', ID);
            s3 = sprintf('GenTx8SpriteAlpha.ppm         !  Sprite File Name\n');

            OverwriteLineInFile(SC_Filename,2,s1);
            OverwriteLineInFile(SC_Filename,3,s2);
            OverwriteLineInFile(SC_Filename,4,s3);
        end
    end   
end

fclose(fileID_input);