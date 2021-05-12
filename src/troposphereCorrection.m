function [troposphereDelay] = troposphereCorrection(tD,receiverPos,satPos)
%TROPOSPHERECORRECTION Summary of this function goes here
%   Detailed explanation goes here

%Define constants
k1=77.604;
k2=382000;
Rd=287.054;
gm=9.784;
g=9.80665;

%Convert the receiver position to geodetic coordinates
receiverPos_llh = ecef2llh(receiverPos);

%Get the satellite position in a ENU reference frame
satPos_enu = ecef2enu(receiverPos,satPos,receiverPos_llh(1),receiverPos_llh(2));

%Get the satellite elevation
[az,el] = enu2AzEl(satPos_enu);

%Compute the obliquity factor
M = 1.001 / sqrt(0.002001 + sind(el)^2);

%Get the averaged meteorological parameters
[parameter] = getParameters(receiverPos_llh(1),tD);
parameter = num2cell(parameter);
[P,T,e,beta,lambda] = deal(parameter{:});

%Compute the vertical delay terms at zero altitude
Tz0_dry = (10^-6 * k1 * Rd * P)/gm;
Tz0_wet = ((10^-6 * k2 * Rd )/((lambda + 1) * gm - (beta * Rd))) * (e/T);


%Compute the vertical delay terms at receiver altitude
H = receiverPos_llh(3);

Tz_dry = (1 - (beta * H)/T)^(g/(Rd*beta)) * Tz0_dry;
Tz_wet = (1 - (beta * H)/T)^(((lambda + 1) * g)/(Rd * beta) - 1) * Tz0_wet;

troposphereDelay = (Tz_dry + Tz_wet) * M;

end

function [parameter] = getParameters(lat,day)

%Define the arrays of the averaged meteorological parameters
avg(1,:) = [1013.25,1017.25,1015.75,1011.75,1013.00]; %P
avg(2,:) = [299.65,294.15,283.15,272.15,263.65]; %T
avg(3,:) = [26.31,21.79,11.66,6.78,4.11]; %e
avg(4,:) = [6.30,6.05,5.58,5.39,4.53] * 0.001; %beta
avg(5,:) = [2.77,3.15,2.57,1.81,1.55]; %lambda

seasonalVar(1,:) = [0.00,-3.75,-2.25,-1.75,-0.50]; %dP
seasonalVar(2,:) = [0.00,7.00,11.00,15.00,14.50]; %dT
seasonalVar(3,:) = [0.00,8.85,7.24,5.36,3.39]; %de
seasonalVar(4,:) = [0.00,0.25,0.32,0.81,0.62] * 0.001; %dbeta
seasonalVar(5,:) = [0.00,0.33,0.46,0.74,0.30]; %dlambda


%Choose the correct values depending on the latitude
if abs(lat) <= 15
    i = 1;
elseif abs(lat) <= 30
    i = 2;
    x1 = 15;
    x2 = 30;
elseif abs(lat) <= 45
    i = 3;
    x1 = 30;
    x2 = 45;
elseif abs(lat) <= 60
    i = 4;
    x1 = 45;
    x2 = 60;
elseif abs(lat) < 75
    i = 5;
    x1 = 60;
    x2 = 75;
else
    i = 6;
end

%Compute the average values
for j = 1:5
    if i == 1
        avgOut(j) = avg(j,1);
        seasonalVarOut(j) = seasonalVar(j,1);
    elseif i == 6
        avgOut(j) = avg(j,5);
        seasonalVarOut(j) = seasonalVar(j,5);
    else
        avgOut(j) = avg(j,i-1) + ((avg(j,i)-avg(j,i-1))/(x2-x1))*(abs(lat)-x1);
        seasonalVarOut(j)=var(j,i-1)+((var(j,i)-var(j,i-1))/(x2-x1))*(abs(lat)-x1);
    end
end

%Compute the meteorological parameters
if lat > 0
    Dmin = 28;
else
    Dmin = 211;
end

for j = 1:5
    parameter(j) = avgOut(j) - seasonalVarOut(j) * cos((2*pi*(day-Dmin))/(365.25));
end

end