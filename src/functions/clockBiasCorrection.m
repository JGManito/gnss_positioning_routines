function [ dtSV ] = clockBiasCorrection(t,af,ecc,sqrtA,E,toc,TGD,prn,debugLevel,fp)
%CLOCKBIASCORRECTION This function computes the clock bias of the SV clock
%   This function implements the clock bias correction described in
%   IS-GPS-200H, page 92. 
%   Inputs: tTX - Time of Transmission
%           tSV - SV PRN Code Phase time at time of transmission
%           af - Vector with clock correction parameters [af0,af1,af2]
%           ecc - Orbital Eccentricity
%           sqrtA - Square root of the Semi-Major Axis
%           E - Mean Anomaly
%           toc - GPS Clock Time
%   Output: dtSV - SV PRN code phase offset, aka clock bias
%

%Define constant F and constants af0 to af2
F = -4.442807633*10^(-10);
af0 = af(1);
af1 = af(2);
af2 = af(3);

%Compute the relativistic correction
dtr = F * ecc * (sqrtA) * sin(E);

%Compute the time difference and correct for the week transition if needed
deltat = t - toc;

if deltat > 302400
    deltat = deltat - 604800;
elseif deltat < -302400
    deltat = deltat + 604800;
end


%Compute the clock bias
dtSV = af0 + af1 * deltat + af2 * deltat^2 + dtr - TGD;

if debugLevel == 1
        fprintf(fp,"PRN %2d: af0=%f, af1=%f, af2=%f, dtr=%+e, deltat=%+e, TGD=%+e, dtSV=%+e\n",prn,af0,af1,af2,dtr,deltat,TGD,dtSV);
end

end

