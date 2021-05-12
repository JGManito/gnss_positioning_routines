function [navMessageFiltered,observationFiltered] = coupleNavObs(t,WN_LSF,observation,navMessage)
%COUPLENAVOBS This function creates an array of navigation messages and
%observations that are coupled to the observation with the same index i and
%valid for the time of the observation
%   Detailed explanation goes here

%Check if any observation has bad health indicators
for i = 1:size(observation,1)
    %Do nothing because the RINEX file used has no health indicators
    %available
end


%Trim the observation array
j = 1;
for i = 1:size(observation,1)
    if sum(observation(i,:)) ~= 0
        observationFiltered(j,:) = observation(i,:);
        j = j+1;
    end
end


%Match the navigation message to the observations
for i = 1:size(observationFiltered,1)
    prnObs = observationFiltered(i,4);
    for prn = 1:size(navMessage,1)
        if prn == prnObs
            navMessageFiltered(i,:) = navMessage(prn,:);
        end
    end
end


%Remove any observations from satellites that aren't present in the
%Navigation Message
temp = navMessageFiltered;
navMessageFiltered = [];
j = 1;
for i = 1:size(temp,1)
    if sum(temp(i,:)) ~= 0
        navMessageFiltered(j,:) = temp(i,:);
        j = j+1;
    end
end

prnNav = navMessageFiltered(:,1);
prnObs = observationFiltered(:,4);

missingIndex = find(ismember(prnObs,prnNav) == 0,1,'first');
observationFiltered(missingIndex,:) = [];



%Check if the navigation message has bad SV health indicator values and
%delete the observation and navigation message in case of bad indicator
i = 1;
n = size(navMessageFiltered,1);
while (i <= n)
    %     disp(i); %DEBUG
    navHealth = navMessageFiltered(i,27);
    if navHealth ~= 0
        navMessageFiltered(i,:) = [];
        observationFiltered(i,:) = [];
    else
        i = i + 1;
    end
    n = size(navMessageFiltered,1);
end


end
