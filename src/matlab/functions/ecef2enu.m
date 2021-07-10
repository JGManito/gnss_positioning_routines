function [ enu ] = ecef2enu(originXYZ,targetXYZ,lat,lon)
%ECEF2ENU Converts coordinates from a cartesian ECEF referential to a local
%ENU referential, with origin in [lat,lon]
%   Detailed explanation goes here


%Get the position difference between target and origin in ECEF referential frame
dxyz = targetXYZ - originXYZ;

%Convert latitude and longitude to radians
lat_r = deg2rad(lat);
lon_r = deg2rad(lon);

%Get the angle for the rotation matrices
alpha = 0.5 * pi - lat_r;
gamma = 0.5 * pi + lon_r;

%Compute the rotation matrices
rotX = [1 0 0;...
        0 cos(alpha) -sin(alpha);...
        0 sin(alpha) cos(alpha)];

rotZ = [cos(gamma) -sin(gamma) 0;...
        sin(gamma) cos(gamma) 0;...
        0 0 1];
    
%Get the coordinates in the ENU reference frame
enu = dxyz * rotZ * rotX;


end

