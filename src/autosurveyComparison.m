clc

%Define reference position value
PositionRefProFlex = [4918524.4824,-791213.3897,3969763.1602];

%Load the ProFlex autosurvey results (obtained by converting the RTKLIB
%results into a CSV file and importing to MATLAB)
load("data\AutoSurvey Results\ProFlex\ProFlex500_03122020_autosurvey.mat");

%The ublox autosurvey record data is obtained by runing the script
%"autosurveyUblox.m"

%% 
for i = 1:size(ProFlex_autosurvey,1)-1
    error_ProFlex_autosurvey(i) = norm(PositionRefProFlex - ProFlex_autosurvey(i,:)); 
end

%%
tPlot = transpose(seconds(0:1:size(ProFlex_autosurvey,1)-1));
figure
hold on
plot(tPlot(1:end-1),error_ubx_svin);
plot(tPlot(1:end-1),error_ProFlex_autosurvey);
xticks(seconds(0:7200:size(error_ubx_svin,2)+1));
xlim([seconds(0) seconds(size(error_ubx_svin,2)+1)]);
xtickangle(45)
xtickformat('hh:mm:ss')
legend('u-blox 6T','ProFlex500')
ylabel('Error (m)')