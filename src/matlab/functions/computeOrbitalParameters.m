function [orbitalParameters] = computeOrbitalParameters(t,navMessage,debugLevel,fp)
%COMPUTEORBITAL This function computes the orbital parameters of a single 
%satellite at the given time t using the ephemeris from the Navigation 
%Message
%   Computes the following parameters (in the order of the output array):
%       1. A - Orbit Semi-Major axis (in meters)
%       2. n - Corrected satellite mean angular velocity
%       3. M - Mean anomaly
%       4. E - Eccentric anomaly
%       5. phi0 - True anomaly
%       6. u - Corrected argument of latitude
%       7. r - Corrected orbital radius
%       8. i - Corrected angle of inclination
%       9. Omega - Corrected longitude of the ascending node


%Define constants
OmegaDot_Earth = 7.2921151467 * 10^-5; %Earth rotation rate as defined in IS-GPS-200K (in radians)
mu_Earth = 3.986005*10^14; %Earth gravitational constant as defined in IS-GPS-200K (in m^3/sec^2)

%Define the precision desired for the eccentric anomaly
prec = 10^-12; %As per RTCM SC 104 (v2.1 p.III-1)

%Split the input variables for easier reading%
Crs = navMessage(9);
deltaN = navMessage(10);
M0 = navMessage(11);
Cuc = navMessage(12);
ecc = navMessage(13);
Cus = navMessage(14);
sqrtA = navMessage(15);
toe = navMessage(16);
Cic = navMessage(17);
Omega0 = navMessage(18);
Cis = navMessage(19);
i0 = navMessage(20);
Crc = navMessage(21);
omega = navMessage(22);
OmegaDot = navMessage(23);
IDOT = navMessage(24);


%Compute the orbit semi-major axis
A = sqrtA^2;

%Compute the corrected satellite mean angular velocity
n = sqrt(mu_Earth/(A^3)) + deltaN;

%Compute the time difference between t and toe
dt = t - toe;

if dt > 302400 %Correct for end of week crossover
    dt = dt - 604800;
elseif dt < -302400
    dt = dt + 604800;
end

%Compute the mean anomaly
M = M0 + n*dt;

%Compute the eccentric anomaly using the iterative method
E = M;
E_ = 0;
while (abs(E - E_) > prec)
    E_ = E;
    E = M + ecc * sin(E); 
end

%Compute the true anomaly
phi0 = atan2(sqrt(1 - ecc^2) * sin(E),cos(E) - ecc);

%Compute the argument of latitude
phi = phi0 + omega;

%Compute the corrected argument of latitude
du = Cuc * cos(2 * phi) + Cus * sin(2 * phi);
u = phi + du;

%Compute the orbital radius
r0 = A * (1 - ecc * cos(E));

%Compute the corrected orbital radius
dr = Crc * cos(2 * phi) + Crs * sin(2 * phi);
r = r0 + dr;

%Compute the corrected angle of inclination
di = Cic * cos(2 * phi) + Cis * sin(2 * phi);
i = i0 + di + IDOT * dt;

%Compute the corrected longitude of the ascending node
Omega = Omega0 + (OmegaDot - OmegaDot_Earth) * dt - OmegaDot_Earth * toe;


%Create the output array
orbitalParameters = [A,n,M,E,phi0,u,r,i,Omega];

if debugLevel == 1
	fprintf(fp,"PRN %2d: ",navMessage(1)); %Move this inside the next function?
    fprintf(fp,"%+f, %+f, %+f, %+f, %+f, %+f, %+f, %+f, %+f, %+f, %+11.6f, %+f, %+f, %+f, %+10.6f\n",A,n,M0,M,E,phi0,phi,du,u,r0,dr,r,di,i,Omega);
end

end

