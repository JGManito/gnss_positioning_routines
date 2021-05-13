clc
clearvars -except gps_observations sbas_observations alpha beta a gps_navigation sbas_navigation
format longg
tic()

%% ---------- Initialization ---------- %%

%Load the dataset
load('MGUE00ARG_R_20202990000_01D_01S_MO.mat');
%Define the reference receiver position
positionRef = [1823327.9138,-4850352.4421,-3709085.5307]; %24h dataset at 1Hz from IGS MGUE marker, Mallarg?e, Argentina
%Set the initial position estimate
initialEstimate_LS = [2756517.813191,-4474878.70094808,-3601428.67912596]; %Buenos Aires, [-34.6,-58.367,0]

%Set the mask angle for 10 degrees
maskAngle = 10;

%Create an array to store the previous observations
observationPrev = zeros(32,2); %Column 1: last pseudorange; Column 2: last phase

%Create an array to store the current navigation message
navMessageCurrent = zeros(32,31);%Create an array to save the most current navigation message for all 32 satellites.
navMessageIndex = 1; %This variable is used as a counter for the current line of the parsed navigation messages

%Set debug level
debugLevel = 1; %Set -1 to disable all output, 0 for results only and 1 for debug



%% ---------- Prepare output file ---------- %%
[fp,errmsg] = fopen('results/pseudorange_partials.txt','w+');


fprintf(fp,"Least-Squares positining using pseudoranges\n\n");

%Output receiver information to file
fprintf(fp,"===== Input Parameters =====\n");
fprintf(fp,"Reference Position: [%f,%f,%f] \n",positionRef(1),positionRef(2),positionRef(3));
fprintf(fp,"Initial Estimate: [%f,%f,%f] \n",initialEstimate_LS(1),initialEstimate_LS(2),initialEstimate_LS(3));
fprintf(fp,"Mask Angle: %.2f degrees \n",maskAngle);


fprintf(fp,"\n\n===== Start of computation =====\n");


%% ---------- Run Least-Squares routine ---------- %%
%Iterate for each observation epoch, t
j = 1; %Counter for the number of iterations
for i=1:16:size(gps_observations,1)
    
    t0 = gps_observations(1,3); %First observation defines the starting time
    
    if i == 1
        dt = gps_observations(17,3) - t0; %Get the time step in seconds
    elseif i <= size(gps_observations,1)-16
        dt = gps_observations(i+16,3) - gps_observations(i,3);
    end
    
    
    t = gps_observations(i,3); %All observations have the same epoch, just use the first one
    WN_LSF = gps_observations(i,2);
    
    %Output time information to file
    fprintf(fp,"Current WN(LSF): %d, Current epoch: %d\n",WN_LSF,t);
    
    
    %Get a slice of the GPS observations array corresponding to the current
    %receiver epoch
    gpsObservationCurrent = gps_observations(i:i+15,:);
    
    %Update the navigation message if there's a new message when
    %the current epoch equals the time of transmission of the navigation
    %message
    [navMessageCurrent,navMessageIndex] = updateNavMessage(gps_navigation,navMessageCurrent,navMessageIndex,t,WN_LSF);
    
    
    %Couple GPS observations with GPS navigation messages
    [navMessageFiltered,observationFiltered] = coupleNavObs(t,WN_LSF,gpsObservationCurrent,navMessageCurrent);
    
    
    %Output satellite information to file
    satList = sprintf('%.0f,',observationFiltered(:,4));
    satList = satList(1:end-1);
    fprintf(fp,"Satellites in view: %s\n",satList);
    
    %Apply the elevation mask
    if t ~= t0
        [observationFiltered,navMessageFiltered] = elevationMask(t,maskAngle,receiverPos_LS(j-1,1:3),observationFiltered,navMessageFiltered,debugLevel,fp);
    end
    
    %Output satellite information to file
    satList = sprintf('%.0f,',navMessageFiltered(:,1));
    satList = satList(1:end-1);
    fprintf(fp,"Satellites in use: %s\n",satList);
    fprintf(fp,"Observations:\n");
    for k=1:size(observationFiltered,1)
        fprintf(fp,"PRN %2d: Pseudorange=%f, CarrierPhase=%f\n",observationFiltered(k,4),observationFiltered(k,5),observationFiltered(k,6));
    end
    
    %Only continue if there's enough satellites
    if size(observationFiltered) >= 4
        
        %----Remove all the modeled errors----%
        [observationFiltered,satPos_tTX,tTX] = removeErrors(observationFiltered,navMessageFiltered,t,alpha,beta,initialEstimate_LS,debugLevel,fp);
        
        %Save the total number and PRN of visible satellites
        nSatsVisible(j) = size(navMessageFiltered,1);
        obsSatUsed(t-t0+1,:) = zeros(1,32);                     %DEBUG
        navSatUsed(t-t0+1,:) = zeros(1,32);                     %DEBUG
        for kk = 1:size(observationFiltered,1)
            currentSVN = observationFiltered(kk,4);
            obsSatUsed(t-t0+1,currentSVN) = currentSVN;
        end
        
        for kk = 1:size(navMessageFiltered,1)                   %DEBUG
            currentSVN = navMessageFiltered(kk,1);              %DEBUG
            navSatUsed(t-t0+1,currentSVN) = currentSVN;         %DEBUG
        end                                                     %DEBUG
        
        if (sum(obsSatUsed(t-t0+1,:)-navSatUsed(t-t0+1,:)) ~=0) %DEBUG
            disp ("nav and obs satellite mismatch!!!");         %DEBUG
        end                                                     %DEBUG
        
        
        %----Least Squares method----%
        [receiverPos_LS(j,:),residuals_LS,GDOP_LS(j),PDOP_LS(j),TDOP_LS(j),H_LS] = leastSquares(observationFiltered,satPos_tTX,initialEstimate_LS,debugLevel,fp);
        
        initialEstimate_LS(1:3) = receiverPos_LS(j,1:3);
        
        positionErrorLS(j) = norm(positionRef - receiverPos_LS(j,1:3));
        errorLS = norm(positionErrorLS(j));
        %fprintf("Position error LS: %f meters\n",errorLS);
        
        if debugLevel >=0
            fprintf(fp,"\nLeast Squares position estimate: [%+15.6f,%+15.6f,%+15.6f]\nReceiver clock error = %+10.6f\n", receiverPos_LS(j,1),receiverPos_LS(j,2),receiverPos_LS(j,3),receiverPos_LS(j,4));
            fprintf(fp,"GDOP = %4.2f, PDOP = %4.2f, TDOP = %4.2f\n", GDOP_LS(j),PDOP_LS(j),TDOP_LS(j));
            fprintf(fp,"Position error: %+15.6f\n",errorLS);
        end
        
        
        j=j+1;
        
        %Separate output information between blocks
        fprintf(fp,"\n--------------------\n");
    end
end

%Close the output file
fclose(fp);