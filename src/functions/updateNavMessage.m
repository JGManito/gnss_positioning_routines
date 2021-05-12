function [navMessageCurrent,navMessageIndex] = updateNavMessage(navMessage,navMessageCurrent,navMessageIndex,tReceiver,WNReceiver)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Get the transmission time and WN of the next navigation message
tTransmission_s = navMessage(navMessageIndex,30);
tTransmission_WN = navMessage(navMessageIndex,3);

%Correct the transmission time for tTX = 999900000 from the CDDIS reference
%dataset
if tTransmission_s == 999900000
    tTransmission_s = 0;
end

%Convert the transmission time and receiver time to seconds since the start
%of GPS time, to account for WN rollover
tTransmission = tTransmission_s + 3600*24*7*tTransmission_WN;
tCurrent = tReceiver + 3600*24*7*WNReceiver;

%Compare the transmission time to the receiver time
while tCurrent >= tTransmission
    %Update the navigation message if there's a new message at current time
    prn = navMessage(navMessageIndex,1);
    navMessageCurrent(prn,:) = navMessage(navMessageIndex,:);
    
    %Increment the index
    if navMessageIndex < size(navMessage,1)
        navMessageIndex = navMessageIndex + 1;
    end
    
    %Check the next transmission to update the cycle condition
    tTransmission_s = navMessage(navMessageIndex,30);
    tTransmission_WN = navMessage(navMessageIndex,3);
    
    %Correct the transmission time for tTX = 999900000 from the CDDIS reference
    %dataset
    if tTransmission_s == 999900000
        tTransmission_s = 0;
    end
    
    tTransmission = tTransmission_s + 3600*24*7*tTransmission_WN;
    
    %Stop when the index has reached the last navigation message
    if navMessageIndex >= size(navMessage,1)
        break;
    end
end


end
