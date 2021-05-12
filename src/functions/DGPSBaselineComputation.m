 function [b,DD,B] = DGPSBaselineComputation(obs1,obs2,satPos,receiverPos1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Save pseudoranges to auxiliary variables
pr1 = obs1(:,5);
pr2 = obs2(:,5);

%Identify the satellite with highest elevation to use as reference for the
%Double Diferences

%Convert the receiver position to geodetic coordinates
receiverPos1_llh = ecef2llh(receiverPos1);

%Get the satellite position in a ENU reference frame
satPos_enu = ecef2enu(receiverPos1,satPos,receiverPos1_llh(1),receiverPos1_llh(2));

%Get the satellite elevation
[~,el] = enu2AzEl(satPos_enu);

%Find the highest elevation satellite
[elMax,n] = max(el);

%Save the observations for that satellite, and remove it from the
%observations
prRefSat1 = pr1(n,:);
prRefSat2 = pr2(n,:);

refSatPos = satPos(n,:);

pr1(n,:)=[];
pr2(n,:)=[];
satPos(n,:) = [];

%Compute the single differences and their direction vectors
SDRef = prRefSat1 - prRefSat2;  %Reference satellite
eRef = (refSatPos - receiverPos1)/norm(refSatPos - receiverPos1);
for i = 1:size(pr1,1)
    SD(i) = pr1(i)-pr2(i);      %Rest of the satellites
    SDe(i,:) = (satPos(i,:) - receiverPos1)/norm(satPos(i,:) - receiverPos1);
end

%Compute the double differences 
for i = 1:size(pr1,1)
    DD(i) = SDRef - SD(i);
    DDe(i,1) = eRef(1) - SDe(i,1);
    DDe(i,2) = eRef(2) - SDe(i,2);
    DDe(i,3) = eRef(3) - SDe(i,3);
end

DD = transpose(DD); %DD is a column vector

%Create the B matrix
B = [DDe(:,1) DDe(:,2) DDe(:,3)];

%Least Squares estimate of the baseline vector
b =inv(B'*B)*B'*DD;

b = transpose(b); %The rest of the code assumes line array for positions

