function [ionosphereDelay] = ionosphereCorrection(positionEstimate,satPosition,t,alpha,beta,prn,debugLevel,fp)
%IONOSPHERECORRECTION This function implements the ionospheric correction
%model as specified in IS-GPS-200K pages 120-122.
%   Inputs:
%       1. Receiver position estimate (LLH)
%       2. Satellite position vector (ECEF, [x,y,z], meters)
%       3. Receiver computed system time (seconds)
%       4. Klobuchar Alpha coefficients (decimal)
%       5. Klobuchar Beta coefficients (decimal)
%
%   Outputs:
%       1. Ionospheric delay (meters)

%Define constants
c = 2.99792458 * 10^8; %WGS-84 definition of c
pi = 3.1415926535898; %Pi value as defined in the IS-GPS-200 specification

%Split input variables for easier reading
alpha0 = alpha(1);
alpha1 = alpha(2);
alpha2 = alpha(3);
alpha3 = alpha(4);

beta0 = beta(1);
beta1 = beta(2);
beta2 = beta(3);
beta3 = beta(4);

%Convert the receiver position to geodetic coordinates
receiverPos_llh = ecef2llh(positionEstimate);

%Compute the satellite Azimuth and Elevation
satPos_ENU = ecef2enu(positionEstimate,satPosition,receiverPos_llh(1),receiverPos_llh(2));
[satAz,satEl] = enu2AzEl(satPos_ENU);
 

%Convert all angular units to semicircles
positionEstimate_sc = rad2semicircles(deg2rad(positionEstimate));
satAz_sc = rad2semicircles(deg2rad(satAz));
satEl_sc = rad2semicircles(deg2rad(satEl));

%-----Klobuchar Model-----%
%Change naming of variables to be consistent with the Klobuchar Model
%variables
lat_u = positionEstimate_sc(1);
lon_u = positionEstimate_sc(2);
gpsTime = t;
A = satAz_sc;
E = satEl_sc;


%Earth's central angle between the user position and the earth projection 
%of ionospheric intersection point
psi = (0.0137 / (E + 0.11)) - 0.022;
    

%Geodetic latitude of the earth projection of the ionospheric intersection point
lat_i = lat_u + psi * cos(A);
if lat_i > 0.416
    lat_i = 0.416;
elseif lat_i < -0.416
    lat_i = -0.416;
end


%Geodetic longitude of the earth projection of the ionospheric intersection point
lon_i = lon_u + (psi * sin(A)) / (cos(lat_i));



%Geomagnetic latitude of the earth projection of the ionospheric intersection point
lat_m = lat_i + 0.064 * cos(lon_i - 1.617);


%Compute local time at the ionospheric intersection point
t = 43200 * lon_i + gpsTime;
if t >= 86400
    t = t - 86400;
elseif t < 0
    t = t + 86400;
end


%Obliquity factor
F = 1 + 16 * (0.53 - E)^3;


%Compute the amplitude of ionospheric delay
AMP = alpha0 + alpha1 * lat_m + alpha2 * lat_m^2 + alpha3 * lat_m^3;
if AMP < 0
    AMP = 0;
end


%Compute the period of ionospheric delay
PER = beta0 + beta1 * lat_m + beta2 * lat_m^2 + beta3 * lat_m^3;
if PER < 72000
    PER = 72000;
end


%Compute the phase of ionospheric delay
X = (2 * pi * (t - 50400))/(PER);


%Compute the ionospheric delay in seconds
if abs(X) < 1.57
    T_iono = F * (5 * 10^-9 + AMP * (1 - (X^2)/2 + (X^4)/24));
else
    T_iono = F * 5 * 10^-9;
end


%Convert the ionospheric delay from seconds to meters
ionosphereDelay = T_iono * c;

if debugLevel == 1
        fprintf(fp,"PRN %2d: A=%f, E=%f, psi=%f, lat_u=%f, lon_u=%f, lat_i=%f, lon_i=%f\n",prn,A,E,psi,lat_u,lon_u,lat_i,lon_i);
        fprintf(fp,"        lat_m=%f, t_IIP=%f, F=%f, AMP=%f, PER=%f, X=%f \n",lat_m,gpsTime,F,AMP,PER,X);
        fprintf(fp,"        IonosphericDelay(s)=%e, IonosphericDelay(m)=%f\n",T_iono,ionosphereDelay);
end

end

