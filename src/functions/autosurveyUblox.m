%---------- UBlox-6T AutoSurvey record processing ----------%

clc
%clear
format longg
tic() %DEBUG
positionRef = [4918525.5233,-791212.0300,3969762.2262];

%% 1 - Load the AutoSurvey record file

%File path of the record
filePath = "data/AutoSurvey Results/COM3_R2_RF2_09112020_SVIN.ubx";

ubx_fp=fopen(filePath,'r','ieee-le');
if (ubx_fp == -1)
    disp("Error reading UBX file")
    return;
end

%% 2 - Parse the AutoSurvey record and extract the SVIN messages

Data = fread(ubx_fp,Inf,'uint8','ieee-le');

[duration,meanX,meanY,meanZ,meanV,nObs,valid,svinActive] = parseSVIN(Data);

%%
meanX = meanX';
meanY = meanY';
meanZ = meanZ';


%% 3 - Compute the error of the AutoSurvey for each observation epoch
for i = 1:size(meanX,1)
    error_ubx_svin(i) = norm(positionRef - [meanX(i),meanY(i),meanZ(i)]);
end

%% 4 - Plot the result
tPlot = transpose(seconds(0:1:size(error_ubx_svin,2)-1));
figure
hold on
plot(tPlot,error_ubx_svin);
xticks(seconds(0:7200:size(error_ubx_svin,2)+1));
xlim([seconds(0) seconds(size(error_ubx_svin,2)+1)]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error (m)')
toc()