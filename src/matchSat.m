function [obs1Out,nav1Out,obs2Out,nav2Out] = matchSat(obs1,nav1,obs2,nav2)
%MATCHSAT This function receives two sets of observations and navigation
%messages and only outputs the data for satellites that are common between
%both sets
%   Detailed explanation goes here

k = 1; %Counter for the size of the output arrays

for i = 1:size(obs1,1)
    for j = 1:size(obs2,1)
        SVN1 = obs1(i,4);
        SVN2 = obs2(j,4);
        if SVN1 == SVN2
            obs1Out(k,:) = obs1(i,:);
            nav1Out(k,:) = nav1(i,:);
            obs2Out(k,:) = obs2(j,:);
            nav2Out(k,:) = nav2(j,:);
            
            k = k + 1;
            break;
        end
    end
end

if size(obs1,1) ~= size (obs1Out,1) || size(obs2,1) ~= size (obs2Out,1)
    disp("Removed satellite!");
    disp("");
end


end

