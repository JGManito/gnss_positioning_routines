function [tTX] = compute_tTX(tRX,pseudorange,clockBias)
%COMPUTE_TTX This function computes the satellite coordinates at
%transmission time and the associated time of transmission in a ECEF
%reference frame


%Define the required constants
c = 2.99792458 * 10^8; %Speed of light as defined in IS-GPS-200K

%Compute the signal flight time
dt = pseudorange/c;

%Compute the time of emission
tTX = tRX - dt - clockBias;