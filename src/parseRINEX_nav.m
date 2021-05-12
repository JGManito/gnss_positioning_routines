function [alpha,beta,a,gps_navigation,sbas_navigation] = parseRINEX_nav(nav_file)
%PARSERINEX_ This function reads the GNSS navigation message file in RINEX format
%and outputs an array with the parsed data for further analysis
%   Detailed explanation goes here


%Header Parsing
%Skip the first lines because they don't provide relevant information.
line = fgetl(nav_file);
while (strcmp(line(1:4),'GPSA') == 0) %Stop skipping when the first ionospheric correction parameters line is found
    line = fgetl(nav_file);
end


%Get the alpha parameters line and separate the values
gpsa = splitNum(line); 
%disp(line); %DEBUG
%disp(gpsa); %DEBUG

alpha0 = exp2float(gpsa(2));
alpha1 = exp2float(gpsa(3));
alpha2 = exp2float(gpsa(4));
alpha3 = exp2float(gpsa(5));

alpha = [alpha0,alpha1,alpha2,alpha3];


%Get the beta parameters line and separate the values
line = fgetl(nav_file);
gpsb = splitNum(line); 
%disp(line); %DEBUG
%disp(gpsb); %DEBUG

beta0 = exp2float(gpsb(2));
beta1 = exp2float(gpsb(3));
beta2 = exp2float(gpsb(4));
beta3 = exp2float(gpsb(5));

beta = [beta0,beta1,beta2,beta3];


%Get the time system corrections line and separate the values
%Currently only GPS is considered. SBAS will be considered after having
%access to real data
line = fgetl(nav_file);
gput = splitNum(line); 
% disp(line); %DEBUG
% disp(gput); %DEBUG

a0 = exp2float(gput(2));
a1 = exp2float(gput(3));
tRef_UTC = str2num(gput{4});
WNRef_UTC = str2num(gput{5});

a=[a0,a1,tRef_UTC,WNRef_UTC]; %System time to UTC correction parameters


%Get the Number of leap seconds since 6-Jan-1980 as transmitted by the GPS almanac
line = fgetl(nav_file);
leapSeconds = split(line); 
% disp(line); %DEBUG
% disp(leapSeconds); %DEBUG

t_LS = leapSeconds{2};


%Skip to the end of the header section
EOH_str1 = '                                                            END OF HEADER';
EOH_str2 = '                                                            END OF HEADER       ';
while (strcmp(line,EOH_str1) ~= 1 && strcmp(line,EOH_str2) ~= 1)
    line = fgetl(nav_file);
%     disp(line) %DEBUG
end

%Parse the rest of the file

i = 1;%Counter for the navigation message index
while(ischar(line))

    line = fgetl(nav_file);
%     disp(line); %DEBUG
    
    if (line(1) == 'G') %Check for the initial line of the navigation message for GPS
        
        %Fix SVN having a space after the G for padding
        if line(2) == ' '
            tempLine = line(1);
            tempLine = strcat(tempLine,line(3:end));
            line = tempLine;
        end
        
        %1st line of the Navigation Message
        navMsg1 = splitNum(line);
        
        
        
        svn_aux = navMsg1{1};
        
        svn = str2num(svn_aux(2:end));
        epochY = str2num(navMsg1{2});
        epochM = str2num(navMsg1{3});
        epochD = str2num(navMsg1{4});
        epochH = str2num(navMsg1{5});
        epochMin = str2num(navMsg1{6});
        epochS = str2num(navMsg1{7});
        
        %Compute GPS ToW and WN
        [epochToW,epochWN,epochWN_LSF] = time2TOW(epochY,epochM,epochD,epochH,epochMin,epochS);
        
        af0 = exp2float(navMsg1{8}); %SV Clock Bias
        af1 = exp2float(navMsg1{9}); %SV Clock drift
        af2 = exp2float(navMsg1{10}); %SV Clock drift rate
        
        af = [af0,af1,af2];
        
        %2nd line of the navigation message
        line = fgetl(nav_file);
%         disp(line); %DEBUG
        navMsg2 = splitNum(line);
        
        IODE = exp2float(navMsg2{1}); %Issue of Data, Ephemeris
        Crs = exp2float(navMsg2{2}); 
        deltaN = exp2float(navMsg2{3});
        Mo = exp2float(navMsg2{4});
        
        %3rd line of the navigation message
        line = fgetl(nav_file);
%         disp(line); %DEBUG
        navMsg3 = splitNum(line);
        
        Cuc = exp2float(navMsg3{1});
        ecc = exp2float(navMsg3{2});
        Cus = exp2float(navMsg3{3});
        sqrtA = exp2float(navMsg3{4});
        
        
        %4th line of the navigation message
        line = fgetl(nav_file);
%         disp(line); %DEBUG
        navMsg4 = splitNum(line);
        
        toe = exp2float(navMsg4{1}); %Time of ephemeris
        Cic = exp2float(navMsg4{2});
        Omega0 = exp2float(navMsg4{3});
        Cis = exp2float(navMsg4{4});
        
        %5th line of the navigation message
        line = fgetl(nav_file);
%         disp(line); %DEBUG
        navMsg5 = splitNum(line);
        
        i0 = exp2float(navMsg5{1});
        Crc = exp2float(navMsg5{2});
        omega = exp2float(navMsg5{3});
        omegaDot = exp2float(navMsg5{4});
        
        %6th line of the navigation message
        line = fgetl(nav_file);
%         disp(line); %DEBUG
        navMsg6 = splitNum(line);
        
        iDot = exp2float(navMsg6{1});
        %L2Codes = exp2float(navMsg6{2}); %Unnecessary?
        WN = exp2float(navMsg6{3});
        %L2Pflag = exp2float(navMsg6{4}); %Unnecessary?
        
        %7th line of the navigation message
        line = fgetl(nav_file);
%         disp(line); %DEBUG
        navMsg7 = splitNum(line);
        
        svAcc = exp2float(navMsg7{1});
        svHealth = exp2float(navMsg7{2});
        TGD = exp2float(navMsg7{3}); %Total Group Dlay
        IODC = exp2float(navMsg7{4});%Issue of Data, Clock
        
        %8th line of the navigation message
        line = fgetl(nav_file);
%         disp(line); %DEBUG
        navMsg8 = splitNum(line);
        
        tTransmission = exp2float(navMsg8{1});     
        fitInterval = 4; %DEBUG
        
        gps_navigation(i,:) = [svn,epochWN,epochWN_LSF,epochToW,af,IODE,Crs,deltaN,Mo,Cuc,ecc,Cus,sqrtA,toe,Cic,Omega0,Cis,i0,...
            Crc,omega,omegaDot,iDot,WN,svAcc,svHealth,TGD,IODC,tTransmission,fitInterval];
        
        i = i+1;
        
    end

sbas_navigation = '';


end


end

