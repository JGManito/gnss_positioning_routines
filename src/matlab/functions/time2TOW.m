function [ToW, WN, WN_LSF] = time2TOW(year,month,day,hour,minute,second)
%TIME2TOW This function converts calendar time to GPS TOW, WN and WN_LSF
%(continuous)
%   Detailed explanation goes here


%Define the start of GPS time as 6/1/1980 at 00:00:00 UTC
GPS_startTime = datetime(1980,1,6,0,0,0,'TimeZone','UTC');


%Convert current time to datetime
GPS_currentTime = datetime(year,month,day,hour,minute,second,'TimeZone','UTC');

%Get Week Number
WN_LSF = split(between(GPS_startTime,GPS_currentTime,'Weeks'),{'Week'});
WN = mod(WN_LSF,1024);

%Get Time of Week
ToW_start = dateshift(GPS_currentTime,'start','week'); %Get the starting time of the current week
dToW = datenum(GPS_currentTime - ToW_start); %Time difference in days
ToW = dToW * 24 * 3600; %Convert to seconds

end

