function [ llh ] = ecef2llh( xyz )
%ECEF2LLH This function converts coordinates in an ECEF reference frame
%(WGS-84 ellipsoid) to Latitude, Longitude and Height using the Heikkinen
%(1982) method
%   Detailed explanation goes here

%Define the WGS-84 ellipsoid constants
a = 6378137;
f = 1/298.257223563;

%Split the input XYZ vector into X, Y and Z for easier reading
x(:) = xyz(:,1);
y(:) = xyz(:,2);
z(:) = xyz(:,3);


%Compute geometrical parameters
b = a * (1 - f);
r = sqrt(x.^2 + y.^2);
e = sqrt(1 - (b^2 / a^2));
ePrime = sqrt((a^2 / b^2) - 1);


%Compute the longitude using Vermeille's (2004) method to avoid
%singularities
if y >= eps
    lon = rad2deg(0.5 * pi - 2 * atan2(x,r+y));
else
    lon = rad2deg(-0.5 * pi + 2 * atan2(x,r-y));
end


%Heikkinen method for latitude and height
F = 54 * b^2 * z.^2;
G = r.^2 + (1 - e^2) .* z.^2 - e^2 * (a^2 - b^2);
c = (e^4 * F .* r.^2) ./ G.^3;
s = nthroot(1 + c + sqrt(c.^2 + 2 * c),3);
P = F ./ (3 * (s + s.^(-1) + 1).^2 .* G.^2);
Q = sqrt(1 + 2 * e^4 * P);
r0 = -(P .* e^2 .* r)/(1 + Q) + sqrt(((a^2) / 2) * (1 + Q.^(-1)) - ...
    ((P * (1 - e^2) .* z.^2) ./ (Q .* (1 + Q))) - ((P .* r.^2)/2));
U = sqrt((r - e^2 * r0).^2 + z.^2);
V = sqrt((r - e^2 * r0).^2 + (1 - e^2)*z.^2);
z0 = (b^2 .* z) / (a * V);

%Finally, calculate height and latitude
h = U .* (1 - (b^2) ./ (a * V));
lat = atan2d(z + ePrime^2 .* z0,r);

%Create the output vector
llh(:,1) = lat;
llh(:,2) = lon;
llh(:,3) = h;

end

