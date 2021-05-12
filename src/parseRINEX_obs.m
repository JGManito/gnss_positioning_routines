function [gps_observations,sbas_observations] = parseRINEX_obs(obs_file)
%PARSERINEX This function reads the GNSS observations file in RINEX format
%and outputs an array with the parsed data for further analysis
%   This function receives the file handle of a RINEX observation file and
%   outputs two arrays:
%       gps_observations(WN,ToW,prn,pseudorange,phase,LLI_flag,doppler,snr)
%       sbas_observations(WN,ToW,prn,pseudorange,phase,LLI_flag,doppler,snr)
%
%       These arrays contain respectively 16 and 8 observations each, with
%       non-existing observations padded with zeros.




line=fgetl(obs_file);
%disp(line) %DEBUG

%Header Parsing
%Skips the header, because it contains no relevant information for the
%current usage
EOH_str1 = '                                                            END OF HEADER';
EOH_str2 = '                                                            END OF HEADER      ';
while (strcmp(line(1:73),EOH_str1) ~= 1 && strcmp(line,EOH_str2) ~= 1)
    line=fgetl(obs_file);
%    disp(line) %DEBUG
end


sbas_observations = 0;

%Saves the observations into an array
nObs = 0; %Counts the number of observations
line=fgetl(obs_file); %Start the cycle
while(ischar(line))
    
    %fprintf("%s\n",line); %DEBUG
    nObs = nObs + 1; %Increment the observation number counter
    
    %Allocate space for the observation arrays.
    gps_observations((nObs-1)*16+1:(nObs-1)*16+16,1:9) = zeros(16,9); %Maximum of 16 possible satellites in GPS constellation
    %sbas_observations((nObs-1)*7+1:(nObs-1)*7+7,1:8) = zeros (7,8); %Currently SBAS consists of 4 WAAS satellites and 3 EGNOS satellites
    
    %Parse the epoch record header
    header=split(line);
    
    %disp(nObs); %DEBUG
    
    tY = str2double(header{2}); %Get observation year
    tM = str2double(header{3}); %Get observation month
    tD = str2double(header{4}); %Get observation day
    th = str2double(header{5}); %Get observation hour
    tm = str2double(header{6}); %Get observation minute
    ts = str2double(header{7}); %Get observation second
    eF = str2double(header{8}); %Get observation epoch flag
    nSats = str2double(header{9}); %Get observation number of satellites
    
    %Compute GPS ToW and WN
    [ToW,WN,WN_LSF] = time2TOW(tY,tM,tD,th,tm,ts);
    %Solve problem with the ToW being parsed as a float
    ToW = str2double(sprintf("%.0f",ToW));
    
    if ToW == 23405
        ToW = 23405;
    end
    
    %Parse the satellite observations
    nSBAS = 1;
    nGPS = 1;
    for i=1:nSats
        line = fgetl(obs_file);
        lengthLine = size(line,2);
        obs = split(line);
        
        %Correct for whitespace in PRN
        if size(obs{1},2) == 1
            obs(1) = strcat(obs(1),'0',obs(2));
            obs(2) = obs(3);
            obs(3) = obs(4);
            obs(4) = obs(5);
            obs(5) = obs(6);
        end
        
        prn = obs{1}; %Get satellite PRN
        
        
        [pseudorange,line] = splitDataFormat(line(4:end),'F14.3');  %Pseudorange
        [~,line] = splitDataFormat(line,'I1');                      %LLI indicator for Pseudorange
        [~,line] = splitDataFormat(line,'I1');                      %Signal Strength
        
        [phase,line] = splitDataFormat(line,'F14.3');               %Phase
        [phaseLLI,line] = splitDataFormat(line,'I1');               %LLI indicator for Pseudorange
        [~,line] = splitDataFormat(line,'I1');                      %Signal Strength
        
        [doppler,line] = splitDataFormat(line,'F14.3');             %Doppler
        [~,line] = splitDataFormat(line,'I1');                      %LLI indicator for Pseudorange
        [~,line] = splitDataFormat(line,'I1');                      %Signal Strength
        
        [snr,~] = splitDataFormat(line,'F14.3');                    %SNR
        if lengthLine >= 69 %If the whitespaces weren't truncated the following variables will exist
            [~,line] = splitDataFormat(line,'I1');                  %LLI indicator for Pseudorange
            [~,line] = splitDataFormat(line,'I1');                  %Signal Strength
        end
        
        
         if (prn(1) == 'S') %Check if SBAS
            %Not implemented
         else %Check if GPS
             prn=str2num(prn(2:3));
             
             gps_observations(nGPS+(nObs-1)*16,:) = [WN,WN_LSF,ToW,prn,pseudorange,phase,phaseLLI,doppler,snr];
             nGPS = nGPS + 1;
         end
        
        end
    if(line == -1)
        pause();
    end
    line=fgetl(obs_file); %Parse new line, here to allow for loop to end on the correct spot
    
end

end

