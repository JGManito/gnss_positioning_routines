function [satXYZ_receptionFrame] = computeSatPosition(orbitalParameters,positionEstimate,prn,debugLevel,fp)
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
satXYZ_emissionFrame = rotZ * rotX * rotZ2 * [r ; 0; 0];



%Convert from a reference system tied to emission time to a reference
%system tied to reception time, common for all measurements

%Calculate the propagation time
dt = norm(satXYZ_emissionFrame - transpose(positionEstimate))/c;
theta = OmegaDot_Earth * dt;

rotZ_Earth = [cos(theta) sin(theta) 0;...
        -sin(theta) cos(theta) 0;...
        0 0 1];

satXYZ_receptionFrame = transpose(rotZ_Earth * satXYZ_emissionFrame);

if debugLevel == 1
   fprintf(fp,"PRN %2d: rotX=[%+f %+f %+f \t rotZ=[%+f %+f %+f \t rotZ2=[%+f %+f %+f \t rotZEarth=[%+f %+f %+f\n",prn,rotX(1,1),rotX(1,2),rotX(1,3),rotZ(1,1),rotZ(1,2),rotZ(1,3),rotZ2(1,1),rotZ2(1,2),rotZ2(1,3),rotZ_Earth(1,1),rotZ_Earth(1,2),rotZ_Earth(1,3)); %Move this inside the next function?
   fprintf(fp,"              %+f %+f %+f \t       %+f %+f %+f \t        %+f %+f %+f \t            %+f %+f %+f\n",rotX(2,1),rotX(2,2),rotX(2,3),rotZ(2,1),rotZ(2,2),rotZ(2,3),rotZ2(2,1),rotZ2(2,2),rotZ2(2,3),rotZ_Earth(2,1),rotZ_Earth(2,2),rotZ_Earth(2,3));
   fprintf(fp,"              %+f %+f %+f]\t       %+f %+f %+f]\t        %+f %+f %+f]\t            %+f %+f %+f]\n",rotX(3,1),rotX(3,2),rotX(3,3),rotZ(3,1),rotZ(3,2),rotZ(3,3),rotZ2(3,1),rotZ2(3,2),rotZ2(3,3),rotZ_Earth(3,1),rotZ_Earth(3,2),rotZ_Earth(3,3));
   fprintf(fp,"        satXYZ_emissionFrame=[%+f,%+f,%+f], dt=%f, theta=%f, satXYZ_receptionFrame=[%+f,%+f,%+f]\n", satXYZ_emissionFrame(1),satXYZ_emissionFrame(2),satXYZ_emissionFrame(3),dt,theta,satXYZ_receptionFrame(1),satXYZ_receptionFrame(2),satXYZ_receptionFrame(3));
end
    

end

