%Data post processing

%% Scatter plot

%Convert ECEF to LLH
receiverPos_LS_LLH = ecef2llh(receiverPos_LS(:,1:3));
receiverPos_WLS_LLH = ecef2llh(receiverPos_WLS(:,1:3));
receiverPos_EKF_LLH = ecef2llh(receiverPos_EKF(:,1:3));
receiverPos_UKF_LLH = ecef2llh(receiverPos_UKF(:,1:3));

%Convert ECEF to ENU
positionRef_LLH=ecef2llh(positionRef);
receiverPos_LS_ENU = ecef2enu(positionRef,receiverPos_LS(:,1:3),positionRef_LLH(1),positionRef_LLH(2));
receiverPos_WLS_ENU = ecef2enu(positionRef,receiverPos_WLS(:,1:3),positionRef_LLH(1),positionRef_LLH(2));
receiverPos_EKF_ENU = ecef2enu(positionRef,receiverPos_EKF(:,1:3),positionRef_LLH(1),positionRef_LLH(2));
receiverPos_UKF_ENU = ecef2enu(positionRef,receiverPos_UKF(:,1:3),positionRef_LLH(1),positionRef_LLH(2));

%2D E-N Scatter Plot
figure
hold on
grid on
scatter(receiverPos_LS_ENU(:,1),receiverPos_LS_ENU(:,2),'+');
scatter(receiverPos_WLS_ENU(:,1),receiverPos_WLS_ENU(:,2),'+');
scatter(receiverPos_EKF_ENU(:,1),receiverPos_EKF_ENU(:,2),'+');
scatter(receiverPos_UKF_ENU(:,1),receiverPos_UKF_ENU(:,2),'+');
maxX = ceil(max(receiverPos_LS_ENU(:,1)))+1;
minX = ceil(min(receiverPos_LS_ENU(:,1)))-1;
maxY = ceil(max(receiverPos_LS_ENU(:,2)))+1;
minY = ceil(min(receiverPos_LS_ENU(:,2)))-1;
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
xlim([minX maxX])
ylim([minY maxY])
xlabel('East (m)')
ylabel('North (m)')
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');


%% Mean Error convergence over time

%The first average error is the error of the first iteration
errorLS_timeDep (1) = norm(positionRef - receiverPos_LS(1,1:3));
errorWLS_timeDep (1) = norm(positionRef - receiverPos_WLS(1,1:3));
errorEKF_timeDep (1) = norm(positionRef - receiverPos_EKF(1,1:3));
errorUKF_timeDep (1) = norm(positionRef - receiverPos_UKF(1,1:3));

parfor i=2:size(receiverPos_LS,1)
    errorLS_timeDep (i) = norm(positionRef - mean(receiverPos_LS(1:i,1:3)));
    errorWLS_timeDep (i) = norm(positionRef - mean(receiverPos_WLS(1:i,1:3)));
    errorEKF_timeDep (i) = norm(positionRef - mean(receiverPos_EKF(1:i,1:3)));
    errorUKF_timeDep (i) = norm(positionRef - mean(receiverPos_UKF(1:i,1:3)));
end

tPlot = transpose(seconds(0:1:size(errorLS_timeDep,2)-1));
figure
hold on
plot(tPlot,errorLS_timeDep,tPlot,errorWLS_timeDep,tPlot,errorEKF_timeDep,tPlot,errorUKF_timeDep);
legend('Least Squares','Weighted Least Squares','Extended Kalman Filter','Unscented Kalman Filter');
xticks(seconds(0:7200:size(errorLS_timeDep,2)));
xlim([seconds(0) seconds(size(errorLS_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error (m)')



%%  LS Error Metrics convergence over time
for i = 1:size(receiverPos_LS,1)
    positionRMS_LS_timeDep(i)   = sqrt(sum((receiverPos_LS(1:i,1)-positionRef(1)).^2 + (receiverPos_LS(1:i,2)-positionRef(2)).^2 + (receiverPos_LS(1:i,3)-positionRef(3)).^2)/i);
    positionDRMS_LS_timeDep(i)  = accMetrics2d('drms',receiverPos_LS(1:i,1:3));
    positionCEP_LS_timeDep(i)   = accMetrics2d('cep',receiverPos_LS(1:i,1:3));
    positionR95_LS_timeDep(i)   = accMetrics2d('r95',receiverPos_LS(1:i,1:3));
    positionMRSE_LS_timeDep(i)  = accMetrics3d('mrse',receiverPos_LS(1:i,1:3));
    positionSEP_LS_timeDep(i)   = accMetrics3d('sep',receiverPos_LS(1:i,1:3));
    positionSAS90_LS_timeDep(i) = accMetrics3d('sas90',receiverPos_LS(1:i,1:3));
end

%%
tPlot = transpose(seconds(0:60:size(positionRMS_LS_timeDep,2)*60));
tPlot = tPlot(1:end-1);
figure
hold on
plot(tPlot,positionDRMS_LS_timeDep,tPlot,positionCEP_LS_timeDep,tPlot,positionR95_LS_timeDep);
legend('DRMS','CEP','R95');
xticks(seconds(0:7200:size(positionRMS_LS_timeDep,2)*60));
xlim([seconds(0) seconds(size(positionRMS_LS_timeDep,2)*60)]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')

figure
hold on
plot(tPlot,positionRMS_LS_timeDep,tPlot,positionMRSE_LS_timeDep,tPlot,positionSEP_LS_timeDep,tPlot,positionSAS90_LS_timeDep);
legend('RMS','MRSE','SEP','SAS90');
xticks(seconds(0:7200:size(positionRMS_LS_timeDep,2)*60));
xlim([seconds(0) seconds(size(positionRMS_LS_timeDep,2)*60)]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')

%%  WLS Error Metrics convergence over time
for i = 1:size(receiverPos_WLS,1)
    positionRMS_WLS_timeDep(i)   = sqrt(sum((receiverPos_WLS(1:i,1)-positionRef(1)).^2 + (receiverPos_WLS(1:i,2)-positionRef(2)).^2 + (receiverPos_WLS(1:i,3)-positionRef(3)).^2)/i);
    positionDRMS_WLS_timeDep(i)  = accMetrics2d('drms',receiverPos_WLS(1:i,1:3));
    positionCEP_WLS_timeDep(i)   = accMetrics2d('cep',receiverPos_WLS(1:i,1:3));
    positionR95_WLS_timeDep(i)   = accMetrics2d('r95',receiverPos_WLS(1:i,1:3));
    positionMRSE_WLS_timeDep(i)  = accMetrics3d('mrse',receiverPos_WLS(1:i,1:3));
    positionSEP_WLS_timeDep(i)   = accMetrics3d('sep',receiverPos_WLS(1:i,1:3));
    positionSAS90_WLS_timeDep(i) = accMetrics3d('sas90',receiverPos_WLS(1:i,1:3));
end

%%
tPlot = transpose(seconds(0:1:size(positionRMS_WLS_timeDep,2)-1));
figure
hold on
plot(tPlot,positionRMS_WLS_timeDep,tPlot,positionDRMS_WLS_timeDep,tPlot,positionCEP_WLS_timeDep,tPlot,positionR95_WLS_timeDep);
legend('RMS','DRMS','CEP','R95');
xticks(seconds(0:7200:size(positionRMS_WLS_timeDep,2)));
xlim([seconds(0) seconds(size(positionRMS_WLS_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')

figure
hold on
plot(tPlot,positionMRSE_WLS_timeDep,tPlot,positionSEP_WLS_timeDep,tPlot,positionSAS90_WLS_timeDep);
legend('MRSE','SEP','SAS90');
xticks(seconds(0:7200:size(positionRMS_WLS_timeDep,2)));
xlim([seconds(0) seconds(size(positionRMS_WLS_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')

%%  LS Error Metrics convergence over time
for i = 1:size(receiverPos_EKF,1)
    positionRMS_EKF_timeDep(i)   = sqrt(sum((receiverPos_EKF(1:i,1)-positionRef(1)).^2 + (receiverPos_EKF(1:i,2)-positionRef(2)).^2 + (receiverPos_EKF(1:i,3)-positionRef(3)).^2)/i);
    positionDRMS_EKF_timeDep(i)  = accMetrics2d('drms',receiverPos_EKF(1:i,1:3));
    positionCEP_EKF_timeDep(i)   = accMetrics2d('cep',receiverPos_EKF(1:i,1:3));
    positionR95_EKF_timeDep(i)   = accMetrics2d('r95',receiverPos_EKF(1:i,1:3));
    positionMRSE_EKF_timeDep(i)  = accMetrics3d('mrse',receiverPos_EKF(1:i,1:3));
    positionSEP_EKF_timeDep(i)   = accMetrics3d('sep',receiverPos_EKF(1:i,1:3));
    positionSAS90_EKF_timeDep(i) = accMetrics3d('sas90',receiverPos_EKF(1:i,1:3));
end

%%
tPlot = transpose(seconds(0:1:size(positionRMS_EKF_timeDep,2)-1));
figure
hold on
plot(tPlot,positionRMS_EKF_timeDep,tPlot,positionDRMS_EKF_timeDep,tPlot,positionCEP_EKF_timeDep,tPlot,positionR95_EKF_timeDep);
legend('RMS','DRMS','CEP','R95');
xticks(seconds(0:7200:size(positionRMS_EKF_timeDep,2)));
xlim([seconds(0) seconds(size(positionRMS_EKF_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')

figure
hold on
plot(tPlot,positionMRSE_EKF_timeDep,tPlot,positionSEP_EKF_timeDep,tPlot,positionSAS90_EKF_timeDep);
legend('MRSE','SEP','SAS90');
xticks(seconds(0:7200:size(positionRMS_EKF_timeDep,2)));
xlim([seconds(0) seconds(size(positionRMS_EKF_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')

%%  UKF Error Metrics convergence over time
for i = 1:size(receiverPos_UKF,1)
    positionRMS_UKF_timeDep(i)   = sqrt(sum((receiverPos_UKF(1:i,1)-positionRef(1)).^2 + (receiverPos_UKF(1:i,2)-positionRef(2)).^2 + (receiverPos_UKF(1:i,3)-positionRef(3)).^2)/i);
    positionDRMS_UKF_timeDep(i)  = accMetrics2d('drms',receiverPos_UKF(1:i,1:3));
    positionCEP_UKF_timeDep(i)   = accMetrics2d('cep',receiverPos_UKF(1:i,1:3));
    positionR95_UKF_timeDep(i)   = accMetrics2d('r95',receiverPos_UKF(1:i,1:3));
    positionMRSE_UKF_timeDep(i)  = accMetrics3d('mrse',receiverPos_UKF(1:i,1:3));
    positionSEP_UKF_timeDep(i)   = accMetrics3d('sep',receiverPos_UKF(1:i,1:3));
    positionSAS90_UKF_timeDep(i) = accMetrics3d('sas90',receiverPos_UKF(1:i,1:3));
end

%%
tPlot = transpose(seconds(0:1:size(positionRMS_UKF_timeDep,2)-1));
figure
hold on
plot(tPlot,positionRMS_UKF_timeDep,tPlot,positionDRMS_UKF_timeDep,tPlot,positionCEP_UKF_timeDep,tPlot,positionR95_UKF_timeDep);
legend('RMS','DRMS','CEP','R95');
xticks(seconds(0:7200:size(positionRMS_UKF_timeDep,2)));
xlim([seconds(0) seconds(size(positionRMS_UKF_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')

figure
hold on
plot(tPlot,positionMRSE_UKF_timeDep,tPlot,positionSEP_UKF_timeDep,tPlot,positionSAS90_UKF_timeDep);
legend('MRSE','SEP','SAS90');
xticks(seconds(0:7200:size(positionRMS_UKF_timeDep,2)));
xlim([seconds(0) seconds(size(positionRMS_UKF_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error Metrics (m)')
%%
tPlot = transpose(seconds(0:1:size(positionErrorLS_time,2)-1));

figure
hold on
plot(tPlot,errorLS_timeDep,tPlot,errorWLS_timeDep,tPlot,errorEKF_timeDep,tPlot,errorUKF_timeDep)
legend('Least-Squares','Weighted Least-Squares','Extended Kalman Filter','Unscented Kalman Filter');
xticks(seconds(0:7200:size(positionErrorLS_time,2)));
xlim([seconds(0) seconds(size(positionErrorLS_time,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Meam error (m)')

%%
for i = 1:size(receiverPos_LS,1)
    fprintf("========== Error Analysis for %d epochs of survey ==========\n",i);

    %-----Least Squares error analysis-----%
    fprintf("----- Least Squares-----\n");

    fprintf("The least-squares error is: %f meters\n",errorLS_timeDep(i*3600));

    positionRMS_LS(i) = sqrt(sum((receiverPos_LS(1:(i*3600),1)-positionRef(1)).^2 + (receiverPos_LS(1:(i*3600),2)-positionRef(2)).^2 + (receiverPos_LS(1:(i*3600),3)-positionRef(3)).^2)/(i*3600));
    fprintf("The least-squares RMS error is: %f meters\n",positionRMS_LS(i));

    positionDRMS_LS(i)  = accMetrics2d('drms',receiverPos_LS(1:(i*3600),1:3));
    positionCEP_LS(i)   = accMetrics2d('cep',receiverPos_LS(1:(i*3600),1:3));
    positionR95_LS(i)   = accMetrics2d('r95',receiverPos_LS(1:(i*3600),1:3));
    positionMRSE_LS(i)  = accMetrics3d('mrse',receiverPos_LS(1:(i*3600),1:3));
    positionSEP_LS(i)   = accMetrics3d('sep',receiverPos_LS(1:(i*3600),1:3));
    positionSAS90_LS(i) = accMetrics3d('sas90',receiverPos_LS(1:(i*3600),1:3));

    fprintf("The least-squares DRMS error is: %f meters\n",positionDRMS_LS(i));
    fprintf("The least-squares CEP error is: %f meters\n",positionCEP_LS(i));
    fprintf("The least-squares R95 error is: %f meters\n",positionR95_LS(i));
    fprintf("The least-squares MRSE error is: %f meters\n",positionMRSE_LS(i));
    fprintf("The least-squares SEP error is: %f meters\n",positionSEP_LS(i));
    fprintf("The least-squares SAS90 error is: %f meters\n\n",positionSAS90_LS(i));


    %-----Weighted Least Squares error analysis-----%
    fprintf("----- Weighted Least Squares-----\n");

    fprintf("The weighted least-squares error is: %f meters\n",errorWLS_timeDep(i*3600));

    positionRMS_WLS(i) = sqrt(sum((receiverPos_WLS(1:(i*3600),1)-positionRef(1)).^2 + (receiverPos_WLS(1:(i*3600),2)-positionRef(2)).^2 + (receiverPos_WLS(1:(i*3600),3)-positionRef(3)).^2)/(i*3600));
    fprintf("The weighted least-squares RMS error is: %f meters\n",positionRMS_WLS(i));

    positionDRMS_WLS(i)  = accMetrics2d('drms',receiverPos_WLS(1:(i*3600),1:3));
    positionCEP_WLS(i)   = accMetrics2d('cep',receiverPos_WLS(1:(i*3600),1:3));
    positionR95_WLS(i)   = accMetrics2d('r95',receiverPos_WLS(1:(i*3600),1:3));
    positionMRSE_WLS(i)  = accMetrics3d('mrse',receiverPos_WLS(1:(i*3600),1:3));
    positionSEP_WLS(i)   = accMetrics3d('sep',receiverPos_WLS(1:(i*3600),1:3));
    positionSAS90_WLS(i) = accMetrics3d('sas90',receiverPos_WLS(1:(i*3600),1:3));

    fprintf("The weighted least-squares DRMS error is: %f meters\n",positionDRMS_WLS(i));
    fprintf("The weighted least-squares CEP error is: %f meters\n",positionCEP_WLS(i));
    fprintf("The weighted least-squares R95 error is: %f meters\n",positionR95_WLS(i));
    fprintf("The weighted least-squares MRSE error is: %f meters\n",positionMRSE_WLS(i));
    fprintf("The weighted least-squares SEP error is: %f meters\n",positionSEP_WLS(i));
    fprintf("The weighted least-squares SAS90 error is: %f meters\n\n",positionSAS90_WLS(i));


    %-----Extended Kalman Filter error analysis-----%
    fprintf("----- Extended Kalman Filter-----\n");
    fprintf("Note: The first %d interations where ignored\n",EKFStart-1);
    fprintf("The Extended Kalman Filter error is: %f meters\n",errorEKF_timeDep(i*3600));

    positionRMS_EKF(i) = sqrt(sum((receiverPos_EKF(EKFStart:(i*3600),1)-positionRef(1)).^2 + (receiverPos_EKF(EKFStart:(i*3600),2)-positionRef(2)).^2 + (receiverPos_EKF(EKFStart:(i*3600),3)-positionRef(3)).^2)/(i*3600));
    fprintf("The Extended Kalman Filter RMS error is: %f meters\n",positionRMS_EKF(i));

    positionDRMS_EKF(i)  = accMetrics2d('drms',receiverPos_EKF(EKFStart:(i*3600),1:3));
    positionCEP_EKF(i)   = accMetrics2d('cep',receiverPos_EKF(EKFStart:(i*3600),1:3));
    positionR95_EKF(i)   = accMetrics2d('r95',receiverPos_EKF(EKFStart:(i*3600),1:3));
    positionMRSE_EKF(i)  = accMetrics3d('mrse',receiverPos_EKF(EKFStart:(i*3600),1:3));
    positionSEP_EKF(i)   = accMetrics3d('sep',receiverPos_EKF(EKFStart:(i*3600),1:3));
    positionSAS90_EKF(i) = accMetrics3d('sas90',receiverPos_EKF(EKFStart:(i*3600),1:3));

    fprintf("The Extended Kalman Filter DRMS error is: %f meters\n",positionDRMS_EKF(i));
    fprintf("The Extended Kalman Filter CEP error is: %f meters\n",positionCEP_EKF(i));
    fprintf("The Extended Kalman Filter R95 error is: %f meters\n",positionR95_EKF(i));
    fprintf("The Extended Kalman Filter MRSE error is: %f meters\n",positionMRSE_EKF(i));
    fprintf("The Extended Kalman Filter SEP error is: %f meters\n",positionSEP_EKF(i));
    fprintf("The Extended Kalman Filter SAS90 error is: %f meters\n\n",positionSAS90_EKF(i));



    %-----Unscented Kalman Filter error analysis-----%
    fprintf("----- Unscented Kalman Filter-----\n");
    fprintf("Note: The first %d interations where ignored\n",UKFStart-1);

    fprintf("The Unscented Kalman Filter error is: %f meters\n",errorUKF_timeDep(i*3600));

    positionRMS_UKF(i) = sqrt(sum((receiverPos_UKF(EKFStart:(i*3600),1)-positionRef(1)).^2 + (receiverPos_UKF(EKFStart:(i*3600),2)-positionRef(2)).^2 + (receiverPos_UKF(EKFStart:(i*3600),3)-positionRef(3)).^2)/(i*3600));
    fprintf("The Unscented Kalman Filter RMS error is: %f meters\n",positionRMS_UKF(i));

    positionDRMS_UKF(i)  = accMetrics2d('drms',receiverPos_UKF(UKFStart:(i*3600),1:3));
    positionCEP_UKF(i)   = accMetrics2d('cep',receiverPos_UKF(UKFStart:(i*3600),1:3));
    positionR95_UKF(i)   = accMetrics2d('r95',receiverPos_UKF(UKFStart:(i*3600),1:3));
    positionMRSE_UKF(i)  = accMetrics3d('mrse',receiverPos_UKF(UKFStart:(i*3600),1:3));
    positionSEP_UKF(i)   = accMetrics3d('sep',receiverPos_UKF(UKFStart:(i*3600),1:3));
    positionSAS90_UKF(i) = accMetrics3d('sas90',receiverPos_UKF(UKFStart:(i*3600),1:3));

    fprintf("The Unscented Kalman Filter DRMS error is: %f meters\n",positionDRMS_UKF(i));
    fprintf("The Unscented Kalman Filter CEP error is: %f meters\n",positionCEP_UKF(i));
    fprintf("The Unscented Kalman Filter R95 error is: %f meters\n",positionR95_UKF(i));
    fprintf("The Unscented Kalman Filter MRSE error is: %f meters\n",positionMRSE_UKF(i));
    fprintf("The Unscented Kalman Filter SEP error is: %f meters\n",positionSEP_UKF(i));
    fprintf("The Unscented Kalman Filter SAS90 error is: %f meters\n\n",positionSAS90_UKF(i));

end

%% Median

%The first average error is the error of the first iteration
errorLS_med_timeDep (1) = norm(positionRef - receiverPos_LS(1,1:3));
errorWLS_med_timeDep (1) = norm(positionRef - receiverPos_WLS(1,1:3));
errorEKF_med_timeDep (1) = norm(positionRef - receiverPos_EKF(1,1:3));
errorUKF_med_timeDep (1) = norm(positionRef - receiverPos_UKF(1,1:3));

parfor i=2:size(receiverPos_LS,1)
    errorLS_med_timeDep (i) = norm(positionRef - median(receiverPos_LS(1:i,1:3)));
    errorWLS_med_timeDep (i) = norm(positionRef - median(receiverPos_WLS(1:i,1:3)));
    errorEKF_med_timeDep (i) = norm(positionRef - median(receiverPos_EKF(1:i,1:3)));
    errorUKF_med_timeDep (i) = norm(positionRef - median(receiverPos_UKF(1:i,1:3)));
end

tPlot = transpose(seconds(0:1:size(errorLS_med_timeDep,2)-1));
figure
hold on
plot(tPlot,errorLS_med_timeDep,tPlot,errorWLS_med_timeDep,tPlot,errorEKF_med_timeDep,tPlot,errorUKF_med_timeDep);
legend('Least Squares','Weighted Least Squares','Extended Kalman Filter','Unscented Kalman Filter');
xticks(seconds(0:7200:size(errorLS_med_timeDep,2)));
xlim([seconds(0) seconds(size(errorLS_med_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error (m)')

%% Weighted mean

     weighted_LS_positions(:,1) = 1./GDOP_LS(:) .* receiverPos_LS(:,1);
     weighted_LS_positions(:,2) = 1./GDOP_LS(:) .* receiverPos_LS(:,2);
     weighted_LS_positions(:,3) = 1./GDOP_LS(:) .* receiverPos_LS(:,3);

     weighted_WLS_positions(:,1) = 1./GDOP_LS(:) .* receiverPos_WLS(:,1);
     weighted_WLS_positions(:,2) = 1./GDOP_LS(:) .* receiverPos_WLS(:,2);
     weighted_WLS_positions(:,3) = 1./GDOP_LS(:) .* receiverPos_WLS(:,3);

     weighted_EKF_positions(:,1) = 1./GDOP_LS(:) .* receiverPos_EKF(:,1);
     weighted_EKF_positions(:,2) = 1./GDOP_LS(:) .* receiverPos_EKF(:,2);
     weighted_EKF_positions(:,3) = 1./GDOP_LS(:) .* receiverPos_EKF(:,3);

     weighted_UKF_positions(:,1) = 1./GDOP_LS(:) .* receiverPos_UKF(:,1);
     weighted_UKF_positions(:,2) = 1./GDOP_LS(:) .* receiverPos_UKF(:,2);
     weighted_UKF_positions(:,3) = 1./GDOP_LS(:) .* receiverPos_UKF(:,3);

for i = 1:size(receiverPos_LS,1)
%      weighted_LS_positions(i,1) = 1/GDOP_LS(i) * receiverPos_LS(i,1);
%      weighted_LS_positions(i,2) = 1/GDOP_LS(i) * receiverPos_LS(i,2);
%      weighted_LS_positions(i,3) = 1/GDOP_LS(i) * receiverPos_LS(i,3);
%
%      weighted_WLS_positions(i,1) = 1/GDOP_LS(i) * receiverPos_WLS(i,1);
%      weighted_WLS_positions(i,2) = 1/GDOP_LS(i) * receiverPos_WLS(i,2);
%      weighted_WLS_positions(i,3) = 1/GDOP_LS(i) * receiverPos_WLS(i,3);
%
%      weighted_EKF_positions(i,1) = 1/GDOP_LS(i) * receiverPos_EKF(i,1);
%      weighted_EKF_positions(i,2) = 1/GDOP_LS(i) * receiverPos_EKF(i,2);
%      weighted_EKF_positions(i,3) = 1/GDOP_LS(i) * receiverPos_EKF(i,3);
%
%      weighted_UKF_positions(i,1) = 1/GDOP_LS(i) * receiverPos_UKF(i,1);
%      weighted_UKF_positions(i,2) = 1/GDOP_LS(i) * receiverPos_UKF(i,2);
%      weighted_UKF_positions(i,3) = 1/GDOP_LS(i) * receiverPos_UKF(i,3);

%     positionWeightedAvgLS(i,:) = [sum(weighted_LS_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_LS_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_LS_positions(1:i,3))/sum(1./GDOP_LS(1:i))];
%     positionWeightedAvgWLS(i,:) = [sum(weighted_WLS_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_WLS_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_WLS_positions(1:i,3))/sum(1./GDOP_LS(1:i))];
%     positionWeightedAvgEKF(i,:) = [sum(weighted_EKF_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_EKF_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_EKF_positions(1:i,3))/sum(1./GDOP_LS(1:i))];
%     positionWeightedAvgUKF(i,:) = [sum(weighted_UKF_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_UKF_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_UKF_positions(1:i,3))/sum(1./GDOP_LS(1:i))];
%
%     errorLS_weightedAvg_timeDep(i) = norm(positionRef - positionWeightedAvgLS(1:i,:));
%     errorWLS_weightedAvg_timeDep(i) = norm(positionRef - positionWeightedAvgWLS(1:i,:));
%     errorEKF_weightedAvg_timeDep(i) = norm(positionRef - positionWeightedAvgEKF(1:i,:));
%     errorUKF_weightedAvg_timeDep(i) = norm(positionRef - positionWeightedAvgUKF(1:i,:));

%Paralel implementation

    errorLS_weightedAvg_timeDep(i) = norm(positionRef - [sum(weighted_LS_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_LS_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_LS_positions(1:i,3))/sum(1./GDOP_LS(1:i))]);
    errorWLS_weightedAvg_timeDep(i) = norm(positionRef - [sum(weighted_WLS_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_WLS_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_WLS_positions(1:i,3))/sum(1./GDOP_LS(1:i))]);
    errorEKF_weightedAvg_timeDep(i) = norm(positionRef - [sum(weighted_EKF_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_EKF_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_EKF_positions(1:i,3))/sum(1./GDOP_LS(1:i))]);
    errorUKF_weightedAvg_timeDep(i) = norm(positionRef - [sum(weighted_UKF_positions(1:i,1))/sum(1./GDOP_LS(1:i)),sum(weighted_UKF_positions(1:i,2))/sum(1./GDOP_LS(1:i)),sum(weighted_UKF_positions(1:i,3))/sum(1./GDOP_LS(1:i))]);

end



tPlot = transpose(seconds(0:1:size(errorLS_med_timeDep,2)-1));
figure
hold on
plot(tPlot,errorLS_weightedAvg_timeDep,tPlot,errorWLS_weightedAvg_timeDep,tPlot,errorEKF_weightedAvg_timeDep,tPlot,errorUKF_weightedAvg_timeDep);
legend('Least Squares','Weighted Least Squares','Extended Kalman Filter','Unscented Kalman Filter');
xticks(seconds(0:7200:size(errorLS_med_timeDep,2)));
xlim([seconds(0) seconds(size(errorLS_med_timeDep,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error (m)')

%% MLE - Normal Distribution

fprintf("\n========== MLE ==========\n");

[mle_LS(1,:)]=mle(receiverPos_LS(:,1));
[mle_LS(2,:)]=mle(receiverPos_LS(:,2));
[mle_LS(3,:)]=mle(receiverPos_LS(:,3));

errorLS_MLE = norm(positionRef - transpose(mle_LS(:,1)));
fprintf("The least-squares error is: %f meters\n",errorLS_MLE);
%fprintf("The least-squares standard deviation is: %f meters\n",errorLS_MLE);

[mle_WLS(1,:)]=mle(receiverPos_WLS(:,1));
[mle_WLS(2,:)]=mle(receiverPos_WLS(:,2));
[mle_WLS(3,:)]=mle(receiverPos_WLS(:,3));

errorWLS_MLE = norm(positionRef - transpose(mle_WLS(:,1)));
fprintf("The least-squares error is: %f meters\n",errorWLS_MLE);
%fprintf("The least-squares standard deviation is: %f meters\n",errorLS_MLE);


[mle_EKF(1,:)]=mle(receiverPos_EKF(:,1));
[mle_EKF(2,:)]=mle(receiverPos_EKF(:,2));
[mle_EKF(3,:)]=mle(receiverPos_EKF(:,3));

errorEKF_MLE = norm(positionRef - transpose(mle_EKF(:,1)));
fprintf("The Extended Kalman Filter error is: %f meters\n",errorEKF_MLE);


[mle_UKF(1,:)]=mle(receiverPos_UKF(:,1));
[mle_UKF(2,:)]=mle(receiverPos_UKF(:,2));
[mle_UKF(3,:)]=mle(receiverPos_UKF(:,3));

errorUKF_MLE = norm(positionRef - transpose(mle_UKF(:,1)));
fprintf("The Unscented Kalman Filter error is: %f meters\n",errorUKF_MLE);


%The first average error is the error of the first iteration
i=1;
[mle_LS(1,1:i)]=receiverPos_LS(1:i,1);
[mle_LS(2,1:i)]=receiverPos_LS(1:i,2);
[mle_LS(3,1:i)]=receiverPos_LS(1:i,3);
errorLS_MLE(i) = norm(positionRef - transpose(mle_LS(:,1)));

[mle_WLS(1,1:i)]=receiverPos_WLS(1:i,1);
[mle_WLS(2,1:i)]=receiverPos_WLS(1:i,2);
[mle_WLS(3,1:i)]=receiverPos_WLS(1:i,3);
errorWLS_MLE(i) = norm(positionRef - transpose(mle_WLS(:,1)));

[mle_EKF(1,1:i)]=receiverPos_EKF(1:i,1);
[mle_EKF(2,1:i)]=receiverPos_EKF(1:i,2);
[mle_EKF(3,1:i)]=receiverPos_EKF(1:i,3);
errorEKF_MLE(i) = norm(positionRef - transpose(mle_EKF(:,1)));

[mle_UKF(1,1:i)]=receiverPos_UKF(1:i,1);
[mle_UKF(2,1:i)]=receiverPos_UKF(1:i,2);
[mle_UKF(3,1:i)]=receiverPos_UKF(1:i,3);
errorUKF_MLE(i) = norm(positionRef - transpose(mle_UKF(:,1)));

for i=2:size(receiverPos_LS,1)
    [mle_LS(1,:),~]=mle(receiverPos_LS(1:i,1));
    [mle_LS(2,:),~]=mle(receiverPos_LS(1:i,2));
    [mle_LS(3,:),~]=mle(receiverPos_LS(1:i,3));

    errorLS_MLE(i) = norm(positionRef - transpose(mle_LS(:,1)));


    [mle_WLS(1,:),~]=mle(receiverPos_WLS(1:i,1));
    [mle_WLS(2,:),~]=mle(receiverPos_WLS(1:i,2));
    [mle_WLS(3,:),~]=mle(receiverPos_WLS(1:i,3));

    errorWLS_MLE(i) = norm(positionRef - transpose(mle_WLS(:,1)));


    [mle_EKF(1,:),~]=mle(receiverPos_EKF(1:i,1));
    [mle_EKF(2,:),~]=mle(receiverPos_EKF(1:i,2));
    [mle_EKF(3,:),~]=mle(receiverPos_EKF(1:i,3));

    errorEKF_MLE(i) = norm(positionRef - transpose(mle_EKF(:,1)));


    [mle_UKF(1,:),~]=mle(receiverPos_UKF(1:i,1));
    [mle_UKF(2,:),~]=mle(receiverPos_UKF(1:i,2));
    [mle_UKF(3,:),~]=mle(receiverPos_UKF(1:i,3));

    errorUKF_MLE(i) = norm(positionRef - transpose(mle_UKF(:,1)));
end
%%
tPlot = transpose(seconds(0:1:size(errorLS_MLE,2)-1));
figure
hold on
plot(tPlot,errorLS_MLE,tPlot,errorWLS_MLE,tPlot,errorEKF_MLE,tPlot,errorUKF_MLE);
legend('Least Squares','Weighted Least Squares','Extended Kalman Filter','Unscented Kalman Filter');
xticks(seconds(0:7200:size(errorLS_MLE,2)));
xlim([seconds(0) seconds(size(errorLS_MLE,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error (m)')


%% MLE - tLocationScale Distribution
%The first average error is the error of the first iteration

for i=1:60 %to solve initial convergence problems
    [mle_LS_tlocscale(1,i)]=receiverPos_LS(i,1);
    [mle_LS_tlocscale(2,i)]=receiverPos_LS(i,2);
    [mle_LS_tlocscale(3,i)]=receiverPos_LS(i,3);
    errorLS_MLE_tlocscale(i) = norm(positionRef - transpose(mle_LS_tlocscale(:,i)));

    [mle_WLS_tlocscale(1,i)]=receiverPos_WLS(i,1);
    [mle_WLS_tlocscale(2,i)]=receiverPos_WLS(i,2);
    [mle_WLS_tlocscale(3,i)]=receiverPos_WLS(i,3);
    errorWLS_MLE_tlocscale(i) = norm(positionRef - transpose(mle_WLS_tlocscale(:,i)));

    [mle_EKF_tlocscale(1,i)]=receiverPos_EKF(i,1);
    [mle_EKF_tlocscale(2,i)]=receiverPos_EKF(i,2);
    [mle_EKF_tlocscale(3,i)]=receiverPos_EKF(i,3);
    errorEKF_MLE_tlocscale(i) = norm(positionRef - transpose(mle_EKF_tlocscale(:,i)));

    [mle_UKF_tlocscale(1,i)]=receiverPos_UKF(i,1);
    [mle_UKF_tlocscale(2,i)]=receiverPos_UKF(i,2);
    [mle_UKF_tlocscale(3,i)]=receiverPos_UKF(i,3);
    errorUKF_MLE_tlocscale(i) = norm(positionRef - transpose(mle_UKF_tlocscale(:,i)));
end

for i=61:size(receiverPos_LS,1)
    TEMP = mle(receiverPos_LS(1:i,1),'distribution','tLocationScale');
    mle_LS_tlocscale(1,:) = TEMP(1);
    TEMP = mle(receiverPos_LS(1:i,2),'distribution','tLocationScale');
    mle_LS_tlocscale(2,:) = TEMP(1);
    TEMP = mle(receiverPos_LS(1:i,3),'distribution','tLocationScale');
    mle_LS_tlocscale(3,:) = TEMP(1);

    errorLS_MLE_tlocscale(i) = norm(positionRef - transpose(mle_LS_tlocscale(:,1)));


    TEMP = mle(receiverPos_WLS(1:i,1),'distribution','tLocationScale');
    mle_WLS_tlocscale(1,:) = TEMP(1);
    TEMP = mle(receiverPos_WLS(1:i,2),'distribution','tLocationScale');
    mle_WLS_tlocscale(2,:) = TEMP(1);
    TEMP = mle(receiverPos_WLS(1:i,3),'distribution','tLocationScale');
    mle_WLS_tlocscale(3,:) = TEMP(1);

    errorWLS_MLE_tlocscale(i) = norm(positionRef - transpose(mle_WLS_tlocscale(:,1)));

    TEMP = mle(receiverPos_EKF(1:i,1),'distribution','tLocationScale');
    mle_EKF_tlocscale(1,:) = TEMP(1);
    TEMP = mle(receiverPos_EKF(1:i,2),'distribution','tLocationScale');
    mle_EKF_tlocscale(2,:) = TEMP(1);
    TEMP = mle(receiverPos_EKF(1:i,3),'distribution','tLocationScale');
    mle_EKF_tlocscale(3,:) = TEMP(1);

    errorEKF_MLE_tlocscale(i) = norm(positionRef - transpose(mle_EKF_tlocscale(:,1)));

    TEMP = mle(receiverPos_UKF(1:i,1),'distribution','tLocationScale');
    mle_UKF_tlocscale(1,:) = TEMP(1);
    TEMP = mle(receiverPos_UKF(1:i,2),'distribution','tLocationScale');
    mle_UKF_tlocscale(2,:) = TEMP(1);
    TEMP = mle(receiverPos_UKF(1:i,3),'distribution','tLocationScale');
    mle_UKF_tlocscale(3,:) = TEMP(1);

    errorUKF_MLE_tlocscale(i) = norm(positionRef - transpose(mle_UKF_tlocscale(:,1)));
end

%% Thresholding

%Sigma-based thresholding
%start at 1h to obtain a modicum of initial data
tic() %DEBUG
n = 2; %n-sigma

std_LS(1) = n*std(receiverPos_LS(1:3600,1));
std_LS(2) = n*std(receiverPos_LS(1:3600,2));
std_LS(3) = n*std(receiverPos_LS(1:3600,3));
mean_LS(1) = mean(receiverPos_LS(1:3600,1));
mean_LS(2) = mean(receiverPos_LS(1:3600,2));
mean_LS(3) = mean(receiverPos_LS(1:3600,3));

receiverPos_LS_threshold = receiverPos_LS(1:3600,1:3);

ii = 3601;
for i=3601:size(receiverPos_LS,1)
    if abs(mean_LS(1) - receiverPos_LS(i,1)) < std_LS(1) && abs(mean_LS(2) - receiverPos_LS(i,2)) < std_LS(2) && abs(mean_LS(3) - receiverPos_LS(i,3)) < std_LS(3)
        receiverPos_LS_threshold(i,:) = receiverPos_LS(i,1:3);
    else
        receiverPos_LS_threshold(i,:) = [NaN,NaN,NaN];
    end

    %Recalculate std and mean
    std_LS(1) = n*std(receiverPos_LS(1:i,1),'omitnan');
    std_LS(2) = n*std(receiverPos_LS(1:i,2),'omitnan');
    std_LS(3) = n*std(receiverPos_LS(1:i,3),'omitnan');
    mean_LS(1) = mean(receiverPos_LS(1:i,1),'omitnan');
    mean_LS(2) = mean(receiverPos_LS(1:i,2),'omitnan');
    mean_LS(3) = mean(receiverPos_LS(1:i,3),'omitnan');

end

std_WLS(1) = n*std(receiverPos_WLS(1:3600,1));
std_WLS(2) = n*std(receiverPos_WLS(1:3600,2));
std_WLS(3) = n*std(receiverPos_WLS(1:3600,3));
mean_WLS(1) = mean(receiverPos_WLS(1:3600,1));
mean_WLS(2) = mean(receiverPos_WLS(1:3600,2));
mean_WLS(3) = mean(receiverPos_WLS(1:3600,3));

receiverPos_WLS_threshold = receiverPos_WLS(1:3600,1:3);

ii = 3601;
for i=3601:size(receiverPos_WLS,1)
    if abs(mean_WLS(1) - receiverPos_WLS(i,1)) < std_WLS(1) && abs(mean_WLS(2) - receiverPos_WLS(i,2)) < std_WLS(2) && abs(mean_WLS(3) - receiverPos_WLS(i,3)) < std_WLS(3)
        receiverPos_WLS_threshold(i,:) = receiverPos_WLS(i,1:3);
    else
        receiverPos_WLS_threshold(i,:) = [NaN,NaN,NaN];
    end

    %Recalculate std and mean
    std_WLS(1) = n*std(receiverPos_WLS(1:i,1),'omitnan');
    std_WLS(2) = n*std(receiverPos_WLS(1:i,2),'omitnan');
    std_WLS(3) = n*std(receiverPos_WLS(1:i,3),'omitnan');
    mean_mean_WWLS(1) = mean(receiverPos_WLS(1:i,1),'omitnan');
    mean_WWLS(2) = mean(receiverPos_WLS(1:i,2),'omitnan');
    mean_WWLS(3) = mean(receiverPos_WLS(1:i,3),'omitnan');

end

std_EKF(1) = n*std(receiverPos_EKF(1:3600,1));
std_EKF(2) = n*std(receiverPos_EKF(1:3600,2));
std_EKF(3) = n*std(receiverPos_EKF(1:3600,3));
mean_EKF(1) = mean(receiverPos_EKF(1:3600,1));
mean_EKF(2) = mean(receiverPos_EKF(1:3600,2));
mean_EKF(3) = mean(receiverPos_EKF(1:3600,3));

receiverPos_EKF_threshold = receiverPos_EKF(1:3600,1:3);

ii = 3601;
for i=3601:size(receiverPos_EKF,1)
    if abs(mean_EKF(1) - receiverPos_EKF(i,1)) < std_EKF(1) && abs(mean_EKF(2) - receiverPos_EKF(i,2)) < std_EKF(2) && abs(mean_EKF(3) - receiverPos_EKF(i,3)) < std_EKF(3)
        receiverPos_EKF_threshold(i,:) = receiverPos_EKF(i,1:3);
    else
        receiverPos_EKF_threshold(i,:) = [NaN,NaN,NaN];
    end

    %Recalculate std and mean
    std_EKF(1) = n*std(receiverPos_EKF(1:i,1),'omitnan');
    std_EKF(2) = n*std(receiverPos_EKF(1:i,2),'omitnan');
    std_EKF(3) = n*std(receiverPos_EKF(1:i,3),'omitnan');
    mean_EKF(1) = mean(receiverPos_EKF(1:i,1),'omitnan');
    mean_EKF(2) = mean(receiverPos_EKF(1:i,2),'omitnan');
    mean_EKF(3) = mean(receiverPos_EKF(1:i,3),'omitnan');

end

std_UKF(1) = n*std(receiverPos_UKF(1:3600,1));
std_UKF(2) = n*std(receiverPos_UKF(1:3600,2));
std_UKF(3) = n*std(receiverPos_UKF(1:3600,3));
mean_UKF(1) = mean(receiverPos_UKF(1:3600,1));
mean_UKF(2) = mean(receiverPos_UKF(1:3600,2));
mean_UKF(3) = mean(receiverPos_UKF(1:3600,3));

receiverPos_UKF_threshold = receiverPos_UKF(1:3600,1:3);

ii = 3601;
for i=3601:size(receiverPos_UKF,1)
    if abs(mean_UKF(1) - receiverPos_UKF(i,1)) < std_UKF(1) && abs(mean_UKF(2) - receiverPos_UKF(i,2)) < std_UKF(2) && abs(mean_UKF(3) - receiverPos_UKF(i,3)) < std_UKF(3)
        receiverPos_UKF_threshold(i,:) = receiverPos_UKF(i,1:3);
    else
        receiverPos_UKF_threshold(i,:) = [NaN,NaN,NaN];
    end

    %Recalculate std and mean
    std_UKF(1) = n*std(receiverPos_UKF(1:i,1),'omitnan');
    std_UKF(2) = n*std(receiverPos_UKF(1:i,2),'omitnan');
    std_UKF(3) = n*std(receiverPos_UKF(1:i,3),'omitnan');
    mean_UKF(1) = mean(receiverPos_UKF(1:i,1),'omitnan');
    mean_UKF(2) = mean(receiverPos_UKF(1:i,2),'omitnan');
    mean_UKF(3) = mean(receiverPos_UKF(1:i,3),'omitnan');

end

%

errorLS_threshold_timeDep(1) = norm(positionRef - receiverPos_LS_threshold(1,:));
errorWLS_threshold_timeDep(1) = norm(positionRef - receiverPos_WLS_threshold(1,:));
errorEKF_threshold_timeDep(1) = norm(positionRef - receiverPos_EKF_threshold(1,:));
errorUKF_threshold_timeDep(1) = norm(positionRef - receiverPos_UKF_threshold(1,:));

for i = 2:size(receiverPos_LS,1)
%     if i == 4001
%         i = 4001;
%     end
    errorLS_threshold_timeDep(i) = norm(positionRef - mean(receiverPos_LS_threshold(1:i,:),'omitnan'));
%     if isnan(errorLS_threshold_timeDep(i))
%         tPlotLS(i) = NaN;
%     end
    
    errorWLS_threshold_timeDep(i) = norm(positionRef - mean(receiverPos_WLS_threshold(1:i,:),'omitnan'));
    errorEKF_threshold_timeDep(i) = norm(positionRef - mean(receiverPos_EKF_threshold(1:i,:),'omitnan'));
    errorUKF_threshold_timeDep(i) = norm(positionRef - mean(receiverPos_UKF_threshold(1:i,:),'omitnan'));
end

%%
tPlot = transpose(seconds(0:1:size(receiverPos_LS_threshold,1)-1));
figure
hold on
plot(tPlot,errorLS_threshold_timeDep,tPlot,errorWLS_threshold_timeDep,tPlot,errorEKF_threshold_timeDep,tPlot,errorUKF_threshold_timeDep);
legend('Least Squares','Weighted Least Squares','Extended Kalman Filter','Unscented Kalman Filter');
xticks(seconds(0:7200:size(receiverPos_LS_threshold,1)));
xlim([seconds(0) seconds(size(receiverPos_LS_threshold,1))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Error (m)')
toc() %DEBUG

%%
%-----------Error analysis after thresholding-----------%
%clc
fprintf("========== Error Analysis after thresholding ==========\n");

%-----Least Squares error analysis-----%
nPointsLS = size(receiverPos_LS_threshold,1) - sum(isnan(receiverPos_LS_threshold(:,1)));
fprintf("----- Least Squares-----\n");
positionAvgLS_threshold(1) = sum(receiverPos_LS_threshold(:,1),'omitnan')/nPointsLS;
positionAvgLS_threshold(2) = sum(receiverPos_LS_threshold(:,2),'omitnan')/nPointsLS;
positionAvgLS_threshold(3) = sum(receiverPos_LS_threshold(:,3),'omitnan')/nPointsLS;
errorLS_threshold = norm(positionRef - positionAvgLS_threshold);
fprintf("The least-squares error is: %f meters\n",errorLS_threshold);

positionRMS_LS_threshold = sqrt(sum((receiverPos_LS_threshold(:,1)-positionRef(1)).^2 + (receiverPos_LS_threshold(:,2)-positionRef(2)).^2 + (receiverPos_LS_threshold(:,3)-positionRef(3)).^2,'omitnan')/nPointsWLS);
fprintf("The least-squares RMS error is: %f meters\n",positionRMS_LS_threshold);

positionDRMS_LS_threshold  = accMetrics2d('drms',receiverPos_LS_threshold(:,1:3));
positionCEP_LS_threshold   = accMetrics2d('cep',receiverPos_LS_threshold(:,1:3));
positionR95_LS_threshold   = accMetrics2d('r95',receiverPos_LS_threshold(:,1:3));
positionMRSE_LS_threshold  = accMetrics3d('mrse',receiverPos_LS_threshold(:,1:3));
positionSEP_LS_threshold   = accMetrics3d('sep',receiverPos_LS_threshold(:,1:3));
positionSAS90_LS_threshold = accMetrics3d('sas90',receiverPos_LS_threshold(:,1:3));

fprintf("The least-squares DRMS error is: %f meters\n",positionDRMS_LS_threshold);
fprintf("The least-squares CEP error is: %f meters\n",positionCEP_LS_threshold);
fprintf("The least-squares R95 error is: %f meters\n",positionR95_LS_threshold);
fprintf("The least-squares MRSE error is: %f meters\n",positionMRSE_LS_threshold);
fprintf("The least-squares SEP error is: %f meters\n",positionSEP_LS_threshold);
fprintf("The least-squares SAS90 error is: %f meters\n\n",positionSAS90_LS_threshold);


%-----Weighted Least Squares error analysis-----%
fprintf("----- Weighted Least Squares-----\n");
nPointsWLS = size(receiverPos_WLS_threshold,1) - sum(isnan(receiverPos_WLS_threshold(:,1)));
positionAvgWLS_threshold(1) = sum(receiverPos_WLS_threshold(:,1),'omitnan')/nPointsWLS;
positionAvgWLS_threshold(2) = sum(receiverPos_WLS_threshold(:,2),'omitnan')/nPointsWLS;
positionAvgWLS_threshold(3) = sum(receiverPos_WLS_threshold(:,3),'omitnan')/nPointsWLS;
errorWLS_threshold = norm(positionRef - positionAvgWLS_threshold);
fprintf("The weighted least-squares error is: %f meters\n",errorWLS_threshold);

positionRMS_WLS_threshold = sqrt(sum((receiverPos_WLS_threshold(:,1)-positionRef(1)).^2 + (receiverPos_WLS_threshold(:,2)-positionRef(2)).^2 + (receiverPos_WLS_threshold(:,3)-positionRef(3)).^2,'omitnan')/nPointsWLS);
fprintf("The weighted least-squares RMS error is: %f meters\n",positionRMS_WLS_threshold);

positionDRMS_WLS_threshold  = accMetrics2d('drms',receiverPos_WLS_threshold(:,1:3));
positionCEP_WLS_threshold   = accMetrics2d('cep',receiverPos_WLS_threshold(:,1:3));
positionR95_WLS_threshold   = accMetrics2d('r95',receiverPos_WLS_threshold(:,1:3));
positionMRSE_WLS_threshold  = accMetrics3d('mrse',receiverPos_WLS_threshold(:,1:3));
positionSEP_WLS_threshold   = accMetrics3d('sep',receiverPos_WLS_threshold(:,1:3));
positionSAS90_WLS_threshold = accMetrics3d('sas90',receiverPos_WLS_threshold(:,1:3));

fprintf("The weighted least-squares DRMS error is: %f meters\n",positionDRMS_WLS_threshold);
fprintf("The weighted least-squares CEP error is: %f meters\n",positionCEP_WLS_threshold);
fprintf("The weighted least-squares R95 error is: %f meters\n",positionR95_WLS_threshold);
fprintf("The weighted least-squares MRSE error is: %f meters\n",positionMRSE_WLS_threshold);
fprintf("The weighted least-squares SEP error is: %f meters\n",positionSEP_WLS_threshold);
fprintf("The weighted least-squares SAS90 error is: %f meters\n\n",positionSAS90_WLS_threshold);


%-----Extended Kalman Filter error analysis-----%
fprintf("----- Extended Kalman Filter-----\n");
fprintf("Note: The first %d interations where ignored\n",EKFStart-1);
nPointsEKF = size(receiverPos_EKF_threshold(EKFStart:end,:),1) - sum(isnan(receiverPos_EKF_threshold(EKFStart:end,1)));
positionAvgEKF_threshold(1) = sum(receiverPos_EKF_threshold(EKFStart:end,1),'omitnan')/nPointsEKF;
positionAvgEKF_threshold(2) = sum(receiverPos_EKF_threshold(EKFStart:end,2),'omitnan')/nPointsEKF;
positionAvgEKF_threshold(3) = sum(receiverPos_EKF_threshold(EKFStart:end,3),'omitnan')/nPointsEKF;
errorEKF_threshold = norm(positionRef - positionAvgEKF_threshold);
fprintf("The Extended Kalman Filter error is: %f meters\n",errorEKF_threshold);

positionRMS_EKF_threshold = sqrt(sum((receiverPos_EKF_threshold(EKFStart:end,1)-positionRef(1)).^2 + (receiverPos_EKF_threshold(EKFStart:end,2)-positionRef(2)).^2 + (receiverPos_EKF_threshold(EKFStart:end,3)-positionRef(3)).^2,'omitnan')/nPointsEKF);
fprintf("The Extended Kalman Filter RMS error is: %f meters\n",positionRMS_EKF_threshold);

positionDRMS_EKF_threshold  = accMetrics2d('drms',receiverPos_EKF_threshold(:,1:3));
positionCEP_EKF_threshold   = accMetrics2d('cep',receiverPos_EKF_threshold(:,1:3));
positionR95_EKF_threshold   = accMetrics2d('r95',receiverPos_EKF_threshold(:,1:3));
positionMRSE_EKF_threshold  = accMetrics3d('mrse',receiverPos_EKF_threshold(:,1:3));
positionSEP_EKF_threshold   = accMetrics3d('sep',receiverPos_EKF_threshold(:,1:3));
positionSAS90_EKF_threshold = accMetrics3d('sas90',receiverPos_EKF_threshold(:,1:3));

fprintf("The Extended Kalman Filter DRMS error is: %f meters\n",positionDRMS_EKF_threshold);
fprintf("The Extended Kalman Filter CEP error is: %f meters\n",positionCEP_EKF_threshold);
fprintf("The Extended Kalman Filter R95 error is: %f meters\n",positionR95_EKF_threshold);
fprintf("The Extended Kalman Filter MRSE error is: %f meters\n",positionMRSE_EKF_threshold);
fprintf("The Extended Kalman Filter SEP error is: %f meters\n",positionSEP_EKF_threshold);
fprintf("The Extended Kalman Filter SAS90 error is: %f meters\n\n",positionSAS90_EKF_threshold);



%-----Unscented Kalman Filter error analysis-----%
fprintf("----- Unscented Kalman Filter-----\n");
fprintf("Note: The first %d interations where ignored\n",UKFStart-1);
nPointsUKF = size(receiverPos_UKF_threshold(UKFStart:end,:),1) - sum(isnan(receiverPos_UKF_threshold(UKFStart:end,1)));
positionAvgUKF_threshold(1) = sum(receiverPos_UKF_threshold(UKFStart:end,1),'omitnan')/nPointsUKF;
positionAvgUKF_threshold(2) = sum(receiverPos_UKF_threshold(UKFStart:end,2),'omitnan')/nPointsUKF;
positionAvgUKF_threshold(3) = sum(receiverPos_UKF_threshold(UKFStart:end,3),'omitnan')/nPointsUKF;
errorUKF_threshold = norm(positionRef - positionAvgUKF_threshold);
fprintf("The Unscented Kalman Filter error is: %f meters\n",errorUKF_threshold);

positionRMS_UKF_threshold = sqrt(sum((receiverPos_UKF_threshold(UKFStart:end,1)-positionRef(1)).^2 + (receiverPos_UKF_threshold(UKFStart:end,2)-positionRef(2)).^2 + (receiverPos_UKF_threshold(UKFStart:end,3)-positionRef(3)).^2,'omitnan')/nPointsUKF);
fprintf("The Unscented Kalman Filter RMS error is: %f meters\n",positionRMS_UKF_threshold);

positionDRMS_UKF_threshold  = accMetrics2d('drms',receiverPos_UKF_threshold(:,1:3));
positionCEP_UKF_threshold   = accMetrics2d('cep',receiverPos_UKF_threshold(:,1:3));
positionR95_UKF_threshold   = accMetrics2d('r95',receiverPos_UKF_threshold(:,1:3));
positionMRSE_UKF_threshold  = accMetrics3d('mrse',receiverPos_UKF_threshold(:,1:3));
positionSEP_UKF_threshold   = accMetrics3d('sep',receiverPos_UKF_threshold(:,1:3));
positionSAS90_UKF_threshold = accMetrics3d('sas90',receiverPos_UKF_threshold(:,1:3));

fprintf("The Unscented Kalman Filter DRMS error is: %f meters\n",positionDRMS_UKF_threshold);
fprintf("The Unscented Kalman Filter CEP error is: %f meters\n",positionCEP_UKF_threshold);
fprintf("The Unscented Kalman Filter R95 error is: %f meters\n",positionR95_UKF_threshold);
fprintf("The Unscented Kalman Filter MRSE error is: %f meters\n",positionMRSE_UKF_threshold);
fprintf("The Unscented Kalman Filter SEP error is: %f meters\n",positionSEP_UKF_threshold);
fprintf("The Unscented Kalman Filter SAS90 error is: %f meters\n\n",positionSAS90_UKF_threshold);

%