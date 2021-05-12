function [metric] = accMetrics2d(metricName,XYZ_pos)
%ACCMETRICS2D This function computes DRMS, 2DRMS, CEP and R95.
%   Detailed explanation goes here

%Remove any NaN values from the arrays
nanIndex = find(isnan(XYZ_pos(:,1)));
XYZ_pos(nanIndex,:) = [];


%The 2D accuracy is measured in the ENU frame. Convert positions to ENU
refPosLLH = ecef2llh(XYZ_pos);
refPos = mean(XYZ_pos);
ENU_pos = ecef2enu(refPos,XYZ_pos,refPosLLH(1),refPosLLH(2));

switch(metricName)
    case 'drms'
        sigmaE = std(ENU_pos(:,1),'omitnan');
        sigmaN = std(ENU_pos(:,2),'omitnan');
        
        metric = sqrt(sigmaE.^2 + sigmaN.^2);
    case '2drms'
        sigmaE = std(ENU_pos(:,1),'omitnan');
        sigmaN = std(ENU_pos(:,2),'omitnan');
        
        metric = 2 * sqrt(sigmaE.^2 + sigmaN.^2);
    case 'cep'       
        nPoints = size(ENU_pos,1);
        E = ENU_pos(:,1);
        N = ENU_pos(:,2);
        
        %Draw a circle of increasing radius, stop when the radius contains
        %50% of the samples
        
        %Set the maximum radius to the most distant measurement
        maxRad = max([max(ENU_pos(:,1)),max(ENU_pos(:,2))]);
        stop = 0;
        
        for r=0.001:0.001:maxRad 
            counter = 0; %Count the number of points that are inside the circle
            for i=1:nPoints
                if sqrt(E(i).^2 + N(i).^2) <= r
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
        
    case 'r95'      
        nPoints = size(ENU_pos,1);
        E = ENU_pos(:,1);
        N = ENU_pos(:,2);
        
        %Draw a circle of increasing radius, stop when the radius contains
        %50% of the samples
        
        %Set the maximum radius to the most distant measurement
        maxRad = max([max(ENU_pos(:,1)),max(ENU_pos(:,2))]);
        stop = 0;
        
        for r=0.001:0.001:maxRad
            counter = 0; %Count the number of points that are inside the circle
            for i=1:nPoints
                if sqrt(E(i).^2 + N(i).^2) <= r
                    counter = counter + 1;
                end
                if counter >= 0.95 * nPoints %If there's 50% of points inside the circle
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

