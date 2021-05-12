function [F,Qk,Xk_est,Pk_est] = unscentedKF(flag,h_0,h_m2,sigmaR,P_initial,observation,satPos_tTX,x0,dt,F,Qk,Xk_est,Pk_est)
%UNSCENTEDKF Summary of this function goes here
%   Detailed explanation goes here

if strcmp('initialize',flag) == 1
    %Add initialization routine
    [F,Qk,Xk_est,Pk_est] = initializeUKF(h_0,h_m2,sigmaR,P_initial,observation,satPos_tTX,x0,dt);
else
    
    %Define constants
    c = 2.99792458 * 10^8; %WGS-84 definition of c
    
    %Split the input arrays for easier reading
    pseudorange = observation(:,5);
    
    %Rename input variables
    Xkprev_p = Xk_est;
    Pkprev_p = Pk_est;
    
    %Compute the observation noise covariance matrix
    nObs = size(observation,1);
    Rk = sigmaR^2 * eye(nObs);
    
    %Factorize the covariance matrix
    sqrtP = chol(Pkprev_p);
    
    %Compute the initial sigma points
    N = size(Xkprev_p,1); %Get the state vector dimension
    tau = 3 - N; %Compute the scaling factor tau
    
    xSigma(:,1) = Xkprev_p;
    for i=1:N
        xSigma(:,i+1) = Xkprev_p + transpose(sqrt(N + tau) * sqrtP(i,:));
        xSigma(:,i+1+N) = Xkprev_p - transpose(sqrt(N + tau) * sqrtP(i,:));
    end
    
    
    %Compute the associated weights
    W(1) = tau / (N + tau);
    for i=1:N
        W(i+1) = 1 / (2 * (N + tau));
        W(i+1+N) = 1 / (2 * (N + tau));
    end
    
    
    %Prediction step
    %Propagate the sigma points through the system dynamic model
    xSigma_m = F * xSigma;
    
    %Compute Xk_m
    Xk_m = zeros(N,1); %Clear the previous value
    for i = 1:(2 * N + 1)
        Xk_m = Xk_m + W(i) * xSigma_m(:,i); %Do the sum
    end
    
    %Compute the error covariance matrix
    Pxk_m = zeros(N,N);
    for i = 1:(2 * N + 1)
        Pxk_m = Pxk_m + W(i) * (xSigma_m(:,i) - Xk_m) * transpose(xSigma_m(:,i) - Xk_m);
    end
    Pxk_m = Pxk_m + Qk;
    
    
    %Update the sigma points after the prediction step
    %Factorize the covariance matrix
    sqrtP = chol(Pxk_m);
    
    %Compute the new sigma points
    N = size(Xk_m,1); %Get the state vector dimension
    tau = 3 - N; %Compute the scaling factor tau
    
    xSigma(:,1) = Xk_m;
    for i=1:N
        xSigma(:,i+1) = Xk_m + transpose(sqrt(N + tau) * sqrtP(i,:));
        xSigma(:,i+1+N) = Xk_m - transpose(sqrt(N + tau) * sqrtP(i,:));
    end
    
    %Filtering step
    %Propagate the sigma points through the observation model
    zSigma_m = h(xSigma,satPos_tTX);
    
    %Compute the average
    M = size(zSigma_m,1);
    Zk_m = zeros(M,1);
    for i = 1:(2 * N + 1)
        Zk_m = Zk_m + W(i) * zSigma_m(:,i); %Do the sum
    end
    
    %Compute the error covariance
    Pzz = zeros(M,M);
    for i = 1:(2 * N + 1)
        Pzz = Pzz +  W(i) * (zSigma_m(:,i) - Zk_m) * transpose(zSigma_m(:,i) - Zk_m);
    end
    Pzz = Pzz + Rk;
    
    %Compute the error cross-variance
    Pxz = zeros(N,M);
    for i = 1:(2 * N + 1)
        Pxz = Pxz + W(i) * (xSigma_m(:,i) - Xk_m) * transpose(zSigma_m(:,i) - Zk_m);
    end
    
    %Compute the Kalman gain
    Kk = Pxz * inv(Pzz);
    
    %Update the state vector and the covariance matrix
    Xk_p = Xk_m + Kk * (pseudorange - Zk_m);
    Pxk_p = Pxk_m - Kk * Pzz * transpose(Kk);
    
    %Rename variables for easier reading
    Xk_est = Xk_p;
    Pk_est = Pxk_p;
    
end
end


function [F,Qk,Xk_est,Pk_est] = initializeUKF(h_0,h_m2,sigmaR,P_initial,observation,satPos_tTX,x0,dt)
%Define constants
c = 2.99792458 * 10^8; %WGS-84 definition of c

%Split the input arrays for easier reading
pseudorange = observation(:,5);


%Compute the initial position estimate
Xkprev_p = [x0(1);x0(2);x0(3);0;0];
%Xkprev_p = [0;0;0;0;0];

%Compute the initial covariance matrix
Pkprev_p = P_initial .* eye(size(Xkprev_p,1));

%Compute the state transition function F
%P model
F = eye(5);
F(4,5) = dt;


%Compute the noise covariance matrix of the dynamics model
q = 0; %Static receiver
q_phi = h_0 * 0.5;
q_f = 2 * pi^2 * h_m2;

a = dt * (c^2 * (q_phi + (q_f * dt^2)/3));
b = dt * (0.5 * (c^2 * q_f * dt));
c = dt * c^2 * q_f;

Qk = [q*dt 0 0 0 0;...
    0 q*dt 0 0 0;...
    0 0 q*dt 0 0;...
    0 0 0 a b;...
    0 0 0 b c];

%Compute the observation noise covariance matrix
nObs = size(observation,1);
Rk = sigmaR^2 * eye(nObs);


%Factorize the covariance matrix
sqrtP = chol(Pkprev_p);


%Compute the initial sigma points
N = size(Xkprev_p,1); %Get the state vector dimension
tau = 3 - N; %Compute the scaling factor tau

xSigma(:,1) = Xkprev_p;
for i=1:N
    xSigma(:,i+1) = Xkprev_p + transpose(sqrt(N + tau) * sqrtP(i,:));
    xSigma(:,i+1+N) = Xkprev_p - transpose(sqrt(N + tau) * sqrtP(i,:));
end


%Compute the associated weights
W(1) = tau / (N + tau);
for i=1:N
    W(i+1) = 1 / (2 * (N + tau));
    W(i+1+N) = 1 / (2 * (N + tau));
end


%Prediction step
%Propagate the sigma points through the system dynamic model
xSigma_m = F * xSigma;

%Compute Xk_m
Xk_m = zeros(N,1); %Clear the previous value
for i = 1:(2 * N + 1)
    Xk_m = Xk_m + W(i) * xSigma_m(:,i); %Do the sum
end

%Compute the error covariance matrix
Pxk_m = zeros(N,N);
for i = 1:(2 * N + 1)
    Pxk_m = Pxk_m + W(i) * (xSigma_m(:,i) - Xk_m) * transpose(xSigma_m(:,i) - Xk_m);
end
Pxk_m = Pxk_m + Qk;


%Update the sigma points after the prediction step
%Factorize the covariance matrix
sqrtP = chol(Pxk_m);

%Compute the new sigma points
N = size(Xk_m,1); %Get the state vector dimension
tau = 3 - N; %Compute the scaling factor tau

xSigma(:,1) = Xk_m;
for i=1:N
    xSigma(:,i+1) = Xk_m + transpose(sqrt(N + tau) * sqrtP(i,:));
    xSigma(:,i+1+N) = Xk_m - transpose(sqrt(N + tau) * sqrtP(i,:));
end

%Filtering step
%Propagate the sigma points through the observation model
zSigma_m = h(xSigma,satPos_tTX);

%Compute the average
M = size(zSigma_m,1);
Zk_m = zeros(M,1);
for i = 1:(2 * N + 1)
    Zk_m = Zk_m + W(i) * zSigma_m(:,i); %Do the sum
end

%Compute the error covariance
Pzz = zeros(M,M);
for i = 1:(2 * N + 1)
    Pzz = Pzz +  W(i) * (zSigma_m(:,i) - Zk_m) * transpose(zSigma_m(:,i) - Zk_m);
end
Pzz = Pzz + Rk;

%Compute the error cross-variance
Pxz = zeros(N,M);
for i = 1:(2 * N + 1)
    Pxz = Pxz + W(i) * (xSigma_m(:,i) - Xk_m) * transpose(zSigma_m(:,i) - Zk_m);
end

%Compute the Kalman gain
Kk = Pxz * inv(Pzz);

%Update the state vector and the covariance matrix
Xk_p = Xk_m + Kk * (pseudorange - Zk_m);
Pxk_p = Pxk_m + Kk * Pzz * transpose(Kk);

%Rename variables for easier reading
Xk_est = Xk_p;
Pk_est = Pxk_p;


end

function [z] = h(x,xSatellite)
%This function computes the observation model of the UKF across all sigma
%points. Note that each of the 2N+1 sigma points contains N rows of data
%corresponding to the state variables.

%Split the input vectors
%P model only
xRec = x(1,:);
yRec = x(2,:);
zRec = x(3,:);
cdtRec = x(4,:);
xSat = xSatellite(:,1);
ySat = xSatellite(:,2);
zSat = xSatellite(:,3);

for j = 1:size(xRec,2)
    for i = 1:size(xSat,1)
        z(i,j) = sqrt((xSat(i) - xRec(j))^2 + (ySat(i) - yRec(j))^2 + (zSat(i) - zRec(j))^2) + cdtRec(j);
    end
end
end
