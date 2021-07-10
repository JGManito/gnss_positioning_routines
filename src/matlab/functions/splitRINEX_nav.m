    function [navMessageSplit] = splitRINEX_nav(t0,WN,navMessage)
%SPLITRINEX_NAV Summary of this function goes here
%   Detailed explanation goes here

fprintf("Nav Message splitter\n");

%Initial values

for prn = 1:32 %For each satellite
    
    match = 0;
    
    for i = 1:size(navMessage,1)
        
        navPRN = navMessage(i,1);
        
        if navPRN == prn %Match the PRN
            navMessageSplit(prn,:) = navMessage(i,:);
            match = 1;
            break;
        end
    end
    
    if match == 0
        navMessageSplit(prn,:) = zeros(1,31);
    end
end

%Observations


%Find the points where to split the navigation message in consistent epoch
%blocks. This will be used to create a time-series-like navigation message

n=1;
start = 0;

for i = 1:size(navMessage,1)
    
    tStep = t0 + n*7200; %2h step
    
    %tEpoch_S = navMessage(i,4);
    tEpoch_S = navMessage(i,30);
    tEpoch_W = navMessage(i,3);
    tEpoch_IODC = navMessage(i,29);
    
    tStepCont = tStep + 3600*24*7*WN;
    tNav = tEpoch_S + 3600*24*7*tEpoch_W;
    
    dt = tStepCont - tNav;
    
    %if dt >=-3600 && dt <=3600 && start == 0 %WORKING!
    if dt >=-3600 && dt <=3600 && start == 0
        iStart(n) = i;
        start = 1;
    end
    
    if dt <-3600 && dt < 0 && start == 1
        iEnd(n) = i-1;
        
        n= n+1;
        iStart(n) = i;
    end
    
    
    navMessageSplit(navPRN+n,:) = navMessage(i,:);
    
end
iEnd(n) = i;

if ~exist('iStart','var')
    aux = zeros(32-size(navMessage,1),31);
    navMessageSplit = [navMessage;aux];
else
    
    %Do the split
    
    for i = 1:n
        navMessageBlock = navMessage(iStart(i):iEnd(i),:);
        navMessageSplit(32*i+1:32*i+32,:) = zeros(32,31);
        for prn = 1:32
            for j = 1:size(navMessageBlock,1)
                navPRN = navMessageBlock(j,1);
                if navPRN == prn
                    navMessageSplit(prn+32*i,:) = navMessageBlock(j,:);
                end
            end
        end
    end
end    
    
    





end

