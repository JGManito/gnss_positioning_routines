function [pseudorangeSmoothed,hatchParameters,observationPrevious] = filterHatch(observation,hatchParameters,observationPrevious,hatchFilterParameter)
%FILTERHATCH This function implements the traditional Hatch filter for
%carrier smoothing of the code pseudorange
%   Detailed explanation goes here

%Split the input arrays for easier reading
SVN = observation(:,4);
pseudorange = observation(:,5);
phase = observation(:,6);


%Get the Hatch Filter weights for the current observations
for i = 1:size(observation)    
    currentSVN = SVN(i);
    W(i,1) = hatchParameters(currentSVN);
    pseudorangeSmoothedPrevious(i,1) = observationPrevious(currentSVN,1);
    phasePrevious(i,1) = observationPrevious(currentSVN,2);
end

%Smooth the pseudoranges using the Hatch Filter, but only if they diverge
%from the pseudorange by more than 50m
pseudorangeSmoothed = W .* pseudorange + (1-W).*(pseudorangeSmoothedPrevious + phase - phasePrevious);

for i=1:size(observation)
    if abs(pseudorangeSmoothed(i) - pseudorange(i)) >= 50
        pseudorangeSmoothed(i) = pseudorange(i);
        currentSVN = SVN(i);
        hatchParameters(currentSVN) = 1; %Reset the Hatch filter parameters if an undetected cycle slip occurs 
    end
end

%Update the filter weight
for i = 1:size(observation)    
    currentSVN = SVN(i);
    if W(i) <= hatchFilterParameter
        hatchParameters(currentSVN) = hatchFilterParameter;
    else
        hatchParameters(currentSVN) = W(i) - hatchFilterParameter;
    end
end


%Update the previous observation array
for i=1:size(observation,1)
    currentSVN = SVN(i);
    observationPrevious(currentSVN,1) = pseudorangeSmoothed(i);
    observationPrevious(currentSVN,2) = phase(i);
end
        

end

