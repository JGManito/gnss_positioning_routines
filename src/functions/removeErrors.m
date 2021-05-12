function [correctedObservation,satPos_tTX,tTX] = removeErrors(observation,navMessage,t,alpha,beta,initialEstimate)
%REMOVEERRORS This function computes all the modeled errors and removes
%them from the observations
%   Detailed explanation goes here

%Define constants
c = 2.99792458 * 10^8; %WGS-84 definition of speed of light
L1_wl = c /(1575.42 * 10^6);

%Define the required precision
prec = 10^-6;

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

%Get the orbital parameters of each satellite
for i = 1:size(observation,1)
    orbitalParameters(i,:) = computeOrbitalParameters(t,navMessage(i,:));
end


%Clock bias error
E = orbitalParameters(:,4);
for i = 1:size(observation,1)
    errClock(i) = clockBiasCorrection(t,af(i,:),ecc(i),sqrtA(i),E(i),toe(i),TGD(i));
    %fprintf("The clock error for PRN%2d is %fs (%fm)\n",prn(i),errClock(i),c*errClock(i));
    %errClock(i) = 0;
end

%Correct the satellite position to obtain the position at time of
%transmission
for i = 1:size(observation,1)
    tTX(i,:) = compute_tTX(t,pseudorange(i),errClock(i));
end


for i = 1:size(observation,1)
    orbitalParameters(i,:) = computeOrbitalParameters(tTX(i),navMessage(i,:));
end

for i = 1:size(observation,1)
    satPos_tTX(i,:) = computeSatPosition(orbitalParameters(i,:),x0);
end


%Tropospheric delay error
for i = 1:size(observation,1)
    [~,~,~,~,~,~,~,tDoY] = TOW2time(t,WN_LSF(i)); %Get the day of the year for the tropospheric correction
    errTropo(i) = troposphereCorrection(tDoY,x0,satPos_tTX(i,:));

end

%Ionospheric delay error
for i = 1:size(observation,1)
    errIono(i) = ionosphereCorrection(x0,satPos_tTX(i,:),tTX,alpha,beta);
end


%Correct the observables
correctedObservation = observation; %Copy all the observations since only pseudorange is corrected in this function
correctedObservation(:,5) = pseudorange + c.*transpose(errClock) - transpose(errTropo) - transpose(errIono);
correctedObservation(:,6) = L1_wl*phase + c.*transpose(errClock) - transpose(errTropo) + transpose(errIono);

end

