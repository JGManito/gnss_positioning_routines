function [correctedObservation,satPos_tTX,tTX] = removeErrors(observation,navMessage,t,alpha,beta,initialEstimate,debugLevel,fp)
%REMOVEERRORS This function computes all the modeled errors and removes
%them from the observations
%   Detailed explanation goes here

%Define constants
c = 2.99792458 * 10^8; %WGS-84 definition of speed of light
L1_wl = c /(1575.42 * 10^6);


%Split the input arrays for easier reading
prn = observation(:,4);
WN_LSF = navMessage(:,3);
toc = navMessage(:,4);
pseudorange = observation(:,5);
phase = observation(:,6);
af = navMessage(:,5:7);
ecc = navMessage(:,13);
sqrtA = navMessage(:,15);
toe = navMessage(:,16);
TGD = navMessage(:,28);

%Define initial position estimate
x0 = [initialEstimate(1),initialEstimate(2),initialEstimate(3)];

if debugLevel == 1
        fprintf(fp,"----- Error correction -----\n");
end

%Get the orbital parameters of each satellite
if debugLevel == 1
        fprintf(fp,"Orbital parameters at reception time for each satellite\n");
        fprintf(fp,"        A                 n          M0         M          E          phi0       phi        du         u          r0                dr           r                 di         i           Omega\n");
end

for i = 1:size(observation,1)
    orbitalParameters(i,:) = computeOrbitalParameters(t,navMessage(i,:),debugLevel,fp);
end


%Clock bias error
E = orbitalParameters(:,4);
if debugLevel == 1
        fprintf(fp,"Clock bias error for each satellite\n");
end
for i = 1:size(observation,1)
    errClock(i) = clockBiasCorrection(t,af(i,:),ecc(i),sqrtA(i),E(i),toe(i),TGD(i),prn(i),debugLevel,fp);
    %fprintf("The clock error for PRN%2d is %fs (%fm)\n",prn(i),errClock(i),c*errClock(i));
    %errClock(i) = 0;
end

%Correct the satellite position to obtain the position at time of
%transmission
if debugLevel == 1
        fprintf(fp,"Time of transmission (in ECEF frame) for each satellite\n");
end
for i = 1:size(observation,1)
    [tTX(i,:),dt(i,:)] = compute_tTX(t,pseudorange(i),errClock(i));
    if debugLevel == 1
        fprintf(fp,"PRN %2d: tTX=%f, tRX=%f, dt=%f\n",prn(i),tTX(i,:),t,dt(i,:));
    end
end

if debugLevel == 1
        fprintf(fp,"Orbital parameters at transmission time for each satellite\n");
        fprintf(fp,"        A                 n          M0         M          E          phi0       phi        du         u          r0                dr           r                 di         i           Omega\n");
end
for i = 1:size(observation,1)
    orbitalParameters(i,:) = computeOrbitalParameters(tTX(i),navMessage(i,:),debugLevel,fp);
end

if debugLevel == 1
        fprintf(fp,"Satellite position at transmission time for each satellite\n");
end
for i = 1:size(observation,1)
    satPos_tTX(i,:) = computeSatPosition(orbitalParameters(i,:),x0,prn(i),debugLevel,fp);
end


%Tropospheric delay error
if debugLevel == 1
        fprintf(fp,"Tropospheric delay error (m) for each satellite (MOPS Model)\n");
end
for i = 1:size(observation,1)
    [~,~,~,~,~,~,~,tDoY] = TOW2time(t,WN_LSF(i)); %Get the day of the year for the tropospheric correction
    errTropo(i) = troposphereCorrection(tDoY,x0,satPos_tTX(i,:),prn(i),debugLevel,fp);
end

%Ionospheric delay error
if debugLevel == 1
        fprintf(fp,"Ionospheric delay error (m) for each satellite (Klobuchar)\n");
end
for i = 1:size(observation,1)
    errIono(i) = ionosphereCorrection(x0,satPos_tTX(i,:),tTX(i,:),alpha,beta,prn(i),debugLevel,fp);
end


%Correct the observables
correctedObservation = observation; %Copy all the observations since only pseudorange is corrected in this function
correctedObservation(:,5) = pseudorange + c.*transpose(errClock) - transpose(errTropo) - transpose(errIono);
correctedObservation(:,6) = L1_wl*phase + c.*transpose(errClock) - transpose(errTropo) + transpose(errIono);

if debugLevel == 1
    fprintf(fp,"Error correction       corrected        observation       phase(m)         clock error(s)  clock error(m)  tropError(m) ionoError(m)\n");
    for i=1:size(observation,1)
        fprintf(fp,"PRN %2d: Pseudorange:   %15.6f, %16.6f, ---------------, %+e, %+15.6f, %+10.6f, %+10.6f\n", prn(i), correctedObservation(i,5), pseudorange(i), transpose(errClock(i)), c*transpose(errClock(i)), transpose(errTropo(i)), transpose(errIono(i)));
        fprintf(fp,"PRN %2d: Carrier Phase: %15.6f, %16.6f, %15.6f, %+e, %+15.6f, %+10.6f, %+10.6f\n", prn(i), correctedObservation(i,6), phase(i), L1_wl*phase(i), transpose(errClock(i)), c*transpose(errClock(i)), transpose(errTropo(i)), transpose(errIono(i)));
    end
end

end

