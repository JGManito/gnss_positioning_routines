function [metric] = accMetrics3d(metricName,XYZ_pos)
%ACCMETRICS3D This function computes MRSE, SEP, SAS90 and SAS99 given the
%error distribution of X, Y and Z coordinates
%   Detailed explanation goes here

%Remove any NaN values from the arrays
nanIndex = find(isnan(XYZ_pos(:,1)));
XYZ_pos(nanIndex,:) = [];

if size(XYZ_pos,1) > 1
    refPos = mean(XYZ_pos);
else
    refPos = XYZ_pos;
end

switch(metricName)
    case 'mrse'
        sigmaX = std(XYZ_pos(:,1));
        sigmaY = std(XYZ_pos(:,2));
        sigmaZ = std(XYZ_pos(:,3));
        
        metric = sqrt(sigmaX.^2 + sigmaY.^2 + sigmaZ.^2);
        
    case 'sep'       
        nPoints = size(XYZ_pos,1);
        dX = XYZ_pos(:,1) - refPos(1);
        dY = XYZ_pos(:,2) - refPos(2);
        dZ = XYZ_pos(:,3) - refPos(3);
        
        %Draw a sphere of increasing radius, stop when the radius contains
        %50% of the samples
        
        %Set the maximum radius to the most distant measurement
        maxRad = max([max(dX),max(dY),max(dZ)]);
        stop = 0;
        
        for r=0.001:0.001:maxRad 
            counter = 0; %Count the number of points that are inside the circle
            for i=1:nPoints
                if sqrt(dX(i).^2 + dY(i).^2 + dZ(i).^2) <= r
                    counter = counter + 1;
                end
                if counter >= 0.5 * nPoints %If there's 50% of points inside the circle
                    stop = 1;
                    break;
                end
            end
            if stop == 1
                break;
            end
        end
        if isempty(r)
            metric = 0;
        else
            metric = r;
        end
        
    case 'sas90'
        nPoints = size(XYZ_pos,1);
        dX = XYZ_pos(:,1) - refPos(1);
        dY = XYZ_pos(:,2) - refPos(2);
        dZ = XYZ_pos(:,3) - refPos(3);
        
        %Draw a sphere of increasing radius, stop when the radius contains
        %90% of the samples
        
        %Set the maximum radius to the most distant measurement
        maxRad = max([max(dX),max(dY),max(dZ)]);
        stop = 0;
        
        for r=0.001:0.001:maxRad 
            counter = 0; %Count the number of points that are inside the circle
            for i=1:nPoints
                if sqrt(dX(i).^2 + dY(i).^2 + dZ(i).^2) <= r
                    counter = counter + 1;
                end
                if counter >= 0.9 * nPoints %If there's 90% of points inside the circle
                    stop = 1;
                    break;
                end
            end
            if stop == 1
                break;
            end
        end
        if isempty(r)
            metric = 0;
        else
            metric = r;
        end
        
    case 'sas99'
        nPoints = size(XYZ_pos,1);
        dX = XYZ_pos(:,1) - refPos(1);
        dY = XYZ_pos(:,2) - refPos(2);
        dZ = XYZ_pos(:,3) - refPos(3);
        
        %Draw a sphere of increasing radius, stop when the radius contains
        %99% of the samples
        
        %Set the maximum radius to the most distant measurement
        maxRad = max([max(dX),max(dY),max(dZ)]);
        stop = 0;
        
        for r=0.001:0.001:maxRad 
            counter = 0; %Count the number of points that are inside the circle
            for i=1:nPoints
                if sqrt(dX(i).^2 + dY(i).^2 + dZ(i).^2) <= r
                    counter = counter + 1;
                end
                if counter >= 0.99 * nPoints %If there's 99% of points inside the circle
                    stop = 1;
                    break;
                end
            end
            if stop == 1
                break;
            end
        end
        if isempty(r)
            metric = 0;
        else
            metric = r;
        end
        
    otherwise
        disp("Warning: unexpected metric requested")
end

end

