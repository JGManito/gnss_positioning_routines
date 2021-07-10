function [ xyz ] = llh2ecef( llh )
%LLH2ECEF Converts geodetic (in degrees) to cartesian coordinates, using WGS-84 as
%the reference ellipsoid
%   Detailed explanation goes here

%Define the WGS-84 ellipsoid constants
a = 6378137;
f = 1/298.257223563;


%Split the input llh vector into lat, lon and h for easier reading
lat = llh(:,1);
lon = llh(:,2);
h = llh(:,3);


%Compute the radius of curvature in the prime vertical
RN = a ./ sqrt(1 - f * (2 - f) * sind(lat).^2);


%Compute the cartesian coordinates
x = (RN + h) .* cosd(lat) .* cosd(lon);
y = (RN + h) .* cosd(lat) .* sind(lon);
z = ((1 - f)^2 .* RN + h) .* sind(lat);

%Create the output vector
xyz (:,1) = x;
xyz (:,2) = y;
xyz (:,3) = z;

end

