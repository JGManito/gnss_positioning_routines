%---------- Differential GPS implementation ----------%
clc
format longg

%% 0 - Set the initial conditions
%Set the initial position estimate
initialEstimate = llh2ecef([38.7166700,-9.1333300,0]); %Lisbon
%rec1Position = [4918524.81,-791212.20,3969762.29]; %Obtida usando PPP
rec1Position = [4918525.55368945, -791211.791994371, 3969762.91803253];% EKF estimate, 24h
rec2Position = [4918532.1188, -791212.5264, 3969754.7230];

%Set the mask angle
maskAngle = 10;

hatchFilterParameter = 0.01;

% % Set the EKF parameters
% h_0 = 2*10^-19; %Temperature-compensated crystal
% h_m2 = 2*10^-20; %Temperature-compensated crystal
% sigmaR = 10;
% P_initial = [10^6,10^6,10^6,9*10^6,9*10^3];

% Set the EKF parameters for the baseline KF
sigmaR_baseline = 10;
%P_initial_baseline = [0.2,0.2,0.2];
P_initial_baseline = [10^2,10^2,10^2];
initial_estimate_baseline = [5,5,5];


%---------- Initialize variables ---------- %

%Create an array to detect cycle slips and to store the Hatch filter
%weights for each satellite
%hatchParameters = [zeros(32,1),ones(32,1)]; %Column 1: last lock; Column 2: Weight
rec1HatchParameters = ones(32,1);
rec2HatchParameters = ones(32,1);

%Create an array to store the previous observations
rec1ObservationPrev = zeros(32,2); %Column 1: last pseudorange; Column 2: last phase
rec2ObservationPrev = zeros(32,2); %Column 1: last pseudorange; Column 2: last phase

rec1NavMessageCurrent = zeros(32,31);%Create an array to save the most current navigation message for all 32 satellites.
rec2NavMessageCurrent = zeros(32,31);%Create an array to save the most current navigation message for all 32 satellites.
rec1NavMessageIndex = 1; %This variable is used as a counter for the current line of the parsed navigation messages
rec2NavMessageIndex = 1; %This variable is used as a counter for the current line of the parsed navigation messages



%% 1 - Load the observations and navigation messages of the two receivers
rec1Data = load('data/7 day dataset/24h Splits/COM3_R2_RF2_25102020.mat');
rec2Data = load('data/7 day dataset/24h Splits/COM4_R1_RF6_25102020.mat');


%Since variables a, alpha and beta are the same for both datasets, just
%differentiate between each receiver's navigation message and observations

alpha = rec1Data.alpha;
beta = rec1Data.beta;

rec1Observations = rec1Data.gps_observations;
rec2Observations = rec2Data.gps_observations;

rec1NavMessage = rec1Data.gps_navigation;
rec2NavMessage = rec2Data.gps_navigation;

clearvars rec1Data rec2Data

%% 2 - Load the receiver1 positions

for i = 1:6
    positionAvgEKF(i,:) = mean(receiverPos_EKF_threshold(1:3600*4*i,:),'omitnan');
    positionAvgUKF(i,:) = mean(receiverPos_UKF_threshold(1:3600*4*i,:),'omitnan');
end



%%
% ---------- Run the simulator ---------- %


for iterIndex = 1:6
    
    %rec1Position = positionAvgEKF(iterIndex,:);
    rec1Position = positionAvgUKF(iterIndex,:);
    
    
    %Iterate for each observation epoch, t
    j = 1; %Counter for the number of iterations
    skipped = 0; %DEBUG
    
    %Copy the initial estimates
    initialEstimate_LS = initialEstimate;
    
    dt = 1;%DEBUG
    
    
    for i=1:16:size(rec1Observations,1)
        t0 = rec1Observations(1,3); %First observation defines the starting time
        
        t = rec1Observations(i,3); %All observations have the same epoch, just use the first one
        WN_LSF = rec1Observations(i,2);
        
        %Get a slice of the GPS observations array corresponding to the current
        %receiver epoch
        rec1ObservationsCurrent = rec1Observations(i:i+15,:);
        rec2ObservationsCurrent = rec2Observations(i:i+15,:);
        
        %Update the navigation message if there's a new message when
        %the current epoch equals the time of transmission of the navigation
        %message
        [rec1NavMessageCurrent,rec1NavMessageIndex] = updateNavMessage(rec1NavMessage,rec1NavMessageCurrent,rec1NavMessageIndex,t,WN_LSF);
        [rec2NavMessageCurrent,rec2NavMessageIndex] = updateNavMessage(rec2NavMessage,rec2NavMessageCurrent,rec2NavMessageIndex,t,WN_LSF);
        
        %Couple GPS observations with GPS navigation messages
        [rec1NavMessageFiltered,rec1ObservationsFiltered] = coupleNavObs(t,WN_LSF,rec1ObservationsCurrent,rec1NavMessageCurrent);
        [rec2NavMessageFiltered,rec2ObservationsFiltered] = coupleNavObs(t,WN_LSF,rec2ObservationsCurrent,rec2NavMessageCurrent);
        
        %Apply the elevation mask
        %if t ~= t0
        [rec1ObservationsFiltered,rec1NavMessageFiltered] = elevationMask(t,maskAngle,initialEstimate,rec1ObservationsFiltered,rec1NavMessageFiltered);
        [rec2ObservationsFiltered,rec2NavMessageFiltered] = elevationMask(t,maskAngle,initialEstimate,rec2ObservationsFiltered,rec2NavMessageFiltered);
        %end
        
        %Filter-out any satellites that arent in use by both receivers
        [rec1ObservationsFiltered,rec1NavMessageFiltered,rec2ObservationsFiltered,rec2NavMessageFiltered] = matchSat(rec1ObservationsFiltered,rec1NavMessageFiltered,rec2ObservationsFiltered,rec2NavMessageFiltered);
        
        
        % Compute the satellite positions and remove the satellite clock bias
        % and TGD error
        %----Remove all the modeled errors----%
        [rec1ObservationsFiltered,satPos_tTX1,tTX1] = removeErrorsDGPS(rec1ObservationsFiltered,rec1NavMessageFiltered,t,alpha,beta,initialEstimate_LS);
        [rec2ObservationsFiltered,satPos_tTX2,tTX2] = removeErrorsDGPS(rec2ObservationsFiltered,rec2NavMessageFiltered,t,alpha,beta,initialEstimate_LS);
        
        %----Implement the Hatch Filter----%
        [rec1HatchParameters] = detectCycleSlip(rec1ObservationsFiltered,rec1HatchParameters);
        [rec2HatchParameters] = detectCycleSlip(rec2ObservationsFiltered,rec2HatchParameters);
        
        [rec1ObservationsFiltered(:,5),rec1HatchParameters,rec1ObservationPrev] = filterHatch(rec1ObservationsFiltered,rec1HatchParameters,rec1ObservationPrev,hatchFilterParameter);
        [rec2ObservationsFiltered(:,5),rec2HatchParameters,rec2ObservationPrev] = filterHatch(rec2ObservationsFiltered,rec2HatchParameters,rec2ObservationPrev,hatchFilterParameter);
        
        
        [b(j,:),DD,B] = DGPSBaselineComputation(rec1ObservationsFiltered,rec2ObservationsFiltered,satPos_tTX1,rec1Position);
        
        %Kalman filtering of the estimate
        if t == t0
            %initialEstimate = b(j,:);
            initialEstimate = initial_estimate_baseline;
            [Pk_prev,Xk_est,Xk_pred] = positionDataKF('initialize',sigmaR_baseline,P_initial_baseline,b(j,:),initialEstimate(1:3));
            [Pk_prev_DD,Xk_est_DD,Xk_pred_DD] = positionDataKF_DD('initialize',sigmaR_baseline,P_initial_baseline,DD,B,initialEstimate(1:3));
        else
            [Pk_prev,Xk_est,Xk_pred] = positionDataKF('recursive',sigmaR_baseline,P_initial_baseline,b(j,:),initialEstimate(1:3),Pk_prev,Xk_est);
            [Pk_prev_DD,Xk_est_DD,Xk_pred_DD] = positionDataKF_DD('recursive',sigmaR_baseline,P_initial_baseline,DD,B,initialEstimate(1:3),Pk_prev_DD,Xk_est_DD);
        end
        
        %For the DD-based KF
        bKF_DD(j,:) = Xk_est_DD;
        rec2EstimateKF_DD(j,:) = rec1Position + bKF_DD(j,:);
        errorEstimateKF_DD(j) = norm(rec2Position - rec2EstimateKF_DD(j,:));
        %distEstimateKF_DD(j) = norm(rec1Position - rec2EstimateKF_DD(j,:));
        
        bKFAvg_DD = mean(bKF_DD(:,:),1);
        rec2EstimateKFAvg_DD(j,:) = rec1Position + bKFAvg_DD;
        errorEstimateKFAvg_DD(j) = norm(rec2Position - rec2EstimateKFAvg_DD(j,:));
        %distEstimateKFAvg_DD(j) = norm(rec1Position - rec2EstimateKFAvg_DD(j,:));
        
        
        
        
        initialEstimate(1:3) = rec2EstimateKF_DD(j,1:3);
        
        
        j=j+1;
        
        
        if t-t0 == 4*3600 %Break after 4h of position determination
            break;
        end
        
    end
    
    %rec2EstimateKF_DD_EKF(iterIndex,:,:) = rec2EstimateKF_DD;
    %rec2EstimateKFAvg_DD_EKF(iterIndex,:,:) = rec2EstimateKFAvg_DD(j-1,:);
    
    rec2EstimateKF_DD_UKF(iterIndex,:,:) = rec2EstimateKF_DD;
    rec2EstimateKFAvg_DD_UKF(iterIndex,:,:) = rec2EstimateKFAvg_DD(j-1,:);
end

%% EKF processing
tPlot = seconds(1:size(rec2EstimateKF_DD_EKF,2));

for i = 1:6
    for j = 1:size(rec2EstimateKF_DD_EKF,2)
        estPosition(i,:) = rec2EstimateKF_DD_EKF(i,j,:);
        errorRec2KF_DD_EKF(i,j) = norm(rec2Position - estPosition);
    end
end
%%
figure
hold on
for i = 1:6
    errorPlot = errorRec2KF_DD_EKF(i,:);
    plot(tPlot,errorPlot)
end
legend('After 4h','After 8h','After 12h','After 16h','After 20h','After 24h')
xticks(seconds(0:1200:size(errorPlot,2)));
xlim([seconds(0) seconds(size(errorPlot,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Absolute error (m)')


%%

fprintf("\n\n\n========== DGPS Results ==========\n");
fprintf("\n");

for i=1:6
    %rec2EstimateAvg(i,:,:) = [mean(rec2EstimateKF_DD_EKF(i,:,1)),mean(rec2EstimateKF_DD_EKF(i,:,2)),mean(rec2EstimateKF_DD_EKF(i,:,3))];
    rec2EstimateAvg(i,:,:) = [mean(rec2EstimateKF_DD_UKF(i,:,1)),mean(rec2EstimateKF_DD_UKF(i,:,2)),mean(rec2EstimateKF_DD_UKF(i,:,3))];
    temp(:) = rec2EstimateAvg(i,:,:);
    rec2EstimateAvgError(i) = norm(rec2Position - temp);

fprintf("Average error after %d hours: %.3f\n",i*4,rec2EstimateAvgError(i));
    
    
end    
    
%% UKF Processing   
tPlot = seconds(1:size(rec2EstimateKF_DD_UKF,2));

for i = 1:6
    for j = 1:size(rec2EstimateKF_DD_UKF,2)
        estPosition(i,:) = rec2EstimateKF_DD_UKF(i,j,:);
        errorRec2KF_DD_UKF(i,j) = norm(rec2Position - estPosition);
    end
end

figure
hold on
for i = 1:6
    errorPlot = errorRec2KF_DD_UKF(i,:);
    plot(tPlot,errorPlot)
end
legend('After 4h','After 8h','After 12h','After 16h','After 20h','After 24h')
xticks(seconds(0:1200:size(errorPlot,2)));
xlim([seconds(0) seconds(size(errorPlot,2))]);
xtickangle(45)
xtickformat('hh:mm:ss')
ylabel('Absolute error (m)')


%%

fprintf("\n\n\n========== DGPS Results ==========\n");
fprintf("\n");
rec2EstimateAvg=[];
rec2EstimateAvgError=[];
temp=[];

for i=1:6
    rec2EstimateAvg(i,:,:) = [mean(rec2EstimateKF_DD_UKF(i,:,1)),mean(rec2EstimateKF_DD_UKF(i,:,2)),mean(rec2EstimateKF_DD_UKF(i,:,3))];
    temp(1,:) = rec2EstimateAvg(i,:,:);
    rec2EstimateAvgError(i) = norm(rec2Position - temp);

fprintf("Average error after %d hours: %.3f\n",i*4,rec2EstimateAvgError(i));
    
    
end        
    
    
%%

for i = 1:6
    errDif(i,:) = errorRec2KF_DD_UKF(i,:) - errorRec2KF_DD_EKF (i,:);
end

plot(tPlot,errDif)
