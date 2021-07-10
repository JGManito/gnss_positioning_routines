function [sortedNav] = sortRINEX_nav(inputNav)
%SORTRINEXX Summary of this function goes here
%   Detailed explanation goes here

%Allocate the output array
sortedNav = zeros(size(inputNav,1),size(inputNav,2));

%Find the initial navigation messages
for i=1:33 %at most 32 satellites in the initial GPS Navigation Message
    epoch = inputNav(i,4);
    
    if epoch > 0 && epoch < 604000 %Find the first navigation message received after the start of data collection
        break;
    end
end

sortedNav(1:i-1,:) = inputNav(1:i-1,:); %Copy the initial message block

sortedNav(i:end,:) = sortrows(inputNav(i:end,:),30);



end

