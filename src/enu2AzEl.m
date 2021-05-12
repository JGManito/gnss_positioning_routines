function [az,el] = enu2AzEl(enu)
%ENU2AZEL Converts local coordinates in a ENU referential to Azimuth and
%Elevation (in degrees)
%   Detailed explanation goes here

%Split the input coordinates for easier reading
e = enu(:,1);
n = enu(:,2);
u = enu(:,3);


%Compute the azimuth
az = atan2d(e,n);

%Compute the elevation
el = atan2d(u,sqrt(n.^2 + e.^2));


end

