function [pseudorange_filtered,navMessage_filtered] = elevationMask(t,maskAngle,receiverPos,pseudorange,navMessage)
%ELEVATIONMASK Function to apply an elevation mask to the observed
%satellites
%   Detailed explanation goes here

k=1;

for i = 1:size(navMessage,1)
    
    %Get the orbital parameters of each satellite
    orbitalParameters(i,:) = computeOrbitalParameters(t,navMessage(i,:));
    
    %Compute the satellite position
    satPos = computeSatPosition(orbitalParameters(i,:),receiverPos);
    
    %Convert the receiver position to geodetic coordinates
    receiverPos_llh = ecef2llh(receiverPos);
    
    %Get the satellite position in a ENU reference frame
    satPos_enu = ecef2enu(receiverPos,satPos,receiverPos_llh(1),receiverPos_llh(2));
    
    %Get the satellite elevation
    [~,el] = enu2AzEl(satPos_enu);
    
    if el > maskAngle
        navMessage_filtered(k,:) = navMessage(i,:);
        pseudorange_filtered(k,:) = pseudorange(i,:);
        k=k+1;
    end
    
    
end

