function [satXYZ] = computeSatPosition(orbitalParameters,positionEstimate)
%COMPUTESATELLITEPOSITION This function computes the satellite position in
%the ECEF referential using it's orbital parameters
%   Detailed explanation goes here

%Define constants
c = 2.99792458 * 10^8; %WGS-84 definition of speed of light
OmegaDot_Earth = 7.2921151467 * 10^-5; %Earth rotation rate as defined in IS-GPS-200K (in radians)

%Extract the required variables from the orbital parameters array
u = orbitalParameters(6);
r = orbitalParameters(7);
i = orbitalParameters(8);
Omega = orbitalParameters(9);


%Calculate the rotation matrices
rotX = [1 0 0;...
        0 cos(-i) sin(-i);...
        0 -sin(-i) cos(-i)];

rotZ = [cos(-Omega) sin(-Omega) 0;...
        -sin(-Omega) cos(-Omega) 0;...
        0 0 1];
    
rotZ2 = [cos(-u) sin(-u) 0;...
        -sin(-u) cos(-u) 0;...
        0 0 1];


%Calculate the satellite position
satXYZ = rotZ * rotX * rotZ2 * [r ; 0; 0];



%Convert from a reference system tied to emission time to a reference
%system tied to reception time, common for all measurements

%Calculate the propagation time
dt = norm(satXYZ - transpose(positionEstimate))/c;
theta = OmegaDot_Earth * dt;

rotZ_Earth = [cos(theta) sin(theta) 0;...
        -sin(theta) cos(theta) 0;...
        0 0 1];

satXYZ = transpose(rotZ_Earth * satXYZ);
    

end

