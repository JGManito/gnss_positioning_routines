function [ hatchParameters ] = detectCycleSlip(observations,hatchParameters)
%DETECTCYCLESLIP This function detects if any cycle slip has occured and
%resets the Hatch Filter weight accordingly.
%   Detailed explanation goes here

for i=1:size(observations,1)
    currentSVN = observations(i,4);
    LLI_flag = observations(i,7);
    
    if LLI_flag ~= 0
        hatchParameters(currentSVN) = 1;
    end
end
        
    
end
