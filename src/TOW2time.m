function [tYear,tMonth,tDay,tWeek,tHour,tMinute,tSecond,tDayofYear] = TOW2time(ToW,WN_LSF)
%TOW2TIME Summary of this function goes here
%   Detailed explanation goes here

%Define the start of GPS time as 6/1/1980 at 00:00:00 UTC
GPS_startTime = datetime(1980,1,6,0,0,0,'TimeZone','UTC');


%Get the current year and week
GPS_time = dateshift(GPS_startTime,'start','week',WN_LSF);

%Get the current day and time of day
GPS_time = dateshift(GPS_time,'start','second',ToW);

%Extract all the time components
tYear = year(GPS_time);
tMonth = month(GPS_time);
tDay = day(GPS_time);
tWeek = week(GPS_time);
tHour = hour(GPS_time);
tMinute = minute(GPS_time);
tSecond = second(GPS_time);

%Get the day of the year 
tDayofYear = day(GPS_time,'dayofyear');

end

