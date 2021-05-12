function [PHI,Qk,Hk,Pk_m,Xk_est,Xk_pred,Pk_pred,Kk,INOV] = extendedKF(flag,h_0,h_m2,sigmaR,P_initial,observation,satPos_tTX,x0,dt,PHI,Qk,Xk_pred,Pk_pred)
%EXTENDEDKF Summary of this function goes here
%   Detailed explanation goes here

if strcmp('initialize',flag) == 1
    %Run the initialization routine
    [PHI,Qk,Hk,Pk_m,Xk_est,Xk_pred,Pk_pred,Kk,INOV] = initializeEKF(h_0,h_m2,sigmaR,P_initial,observation,satPos_tTX,x0,dt);
else
    
    %Define constants
    c = 2.99792458 * 10^8; %WGS-84 definition of c
    
    %Split the input arrays for easier reading
    pseudorange = observation(:,5);

    
    %Update the indexes of the Extended Kalman Filter arrays
    Xk_m = Xk_pred;
    Pk_m = Pk_pred;
    
    %Compute the observation noise covariance matrix
    nObs = size(observation,1);
    Rk = sigmaR^2 * eye(nObs);
    
    %Compute the observation matrix
    for i = 1:nObs
        r = sqrt((satPos_tTX(i,1) - Xk_m(1))^2 + (satPos_tTX(i,2) - Xk_m(2))^2 + (satPos_tTX(i,3) - Xk_m(3))^2);
        Hk(i,1) = -(satPos_tTX(i,1) - Xk_m(1)) / r;
        Hk(i,2) = -(satPos_tTX(i,2) - Xk_m(2)) / r;
        Hk(i,3) = -(satPos_tTX(i,3) - Xk_m(3)) / r;
        Hk(i,4) = 1;
        Hk(i,5) = 0;
    end
    
    %Compute the Kalman gain matrix
    Kk = Pk_m * Hk' * inv((Hk * Pk_m * Hk' + Rk));
    
    %Update the position estimate
    for i = 1:nObs
        h(i) = sqrt((satPos_tTX(i,1) - Xk_m(1))^2 + (satPos_tTX(i,2) - Xk_m(2))^2 + (satPos_tTX(i,3) - Xk_m(3))^2) + Xk_m(4);
    end
    
    INOV = pseudorange - h';
    
    Xk_p(:,1) = Xk_m + Kk * (pseudorange - h');
    
    
    %Update the error covariance matrix
    I = eye(5);
    Pk_p = (I - Kk * Hk) * Pk_m * (I - Kk * Hk)' + Kk * Rk * Kk';
    
    %Predict the state vector
    Xk_pred(:,1) = PHI * Xk_p;
    
    %Predict the error covariance matrix
    Pk_pred(:,:,1) = PHI * Pk_p * PHI' + Qk;
    
    %Rename the state estimate variable for easier reading
    Xk_est = Xk_p;
    
    
    
end
end

function [PHI,Qk,Hk,Pk_m,Xk_est,Xk_pred,Pk_pred,Kk,INOV] = initializeEKF(h_0,h_m2,sigmaR,P_initial,observation,satPos_tTX,x0,dt)

%Define constants
c = 2.99792458 * 10^8; %WGS-84 definition of c

%Split the input arrays for easier reading
pseudorange = observation(:,5);


%Compute the initial position estimate
Xk_m = [x0(1);x0(2);x0(3);0;0];

%Compute the State Transition Matrix
PHI = eye(5);
PHI(4,5) = dt;

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

%Compute the observation matrix
for i = 1:nObs
    r = sqrt((satPos_tTX(i,1) - Xk_m(1))^2 + (satPos_tTX(i,2) - Xk_m(2))^2 + (satPos_tTX(i,3) - Xk_m(3))^2);
    Hk(i,1) = -(satPos_tTX(i,1) - Xk_m(1)) / r;
    Hk(i,2) = -(satPos_tTX(i,2) - Xk_m(2)) / r;
    Hk(i,3) = -(satPos_tTX(i,3) - Xk_m(3)) / r;
    Hk(i,4) = 1;
    Hk(i,5) = 0;
end


%Compute the error covariance matrix
%Pk_m = P_initial * eye(5);
Pk_m = PHI * (P_initial .* eye(5)) * PHI' + Qk;

%Compute the Kalman gain matrix
Kk = Pk_m * Hk' * inv((Hk * Pk_m * Hk' + Rk));

%Update the position estimate
for i = 1:nObs
    h(i) = sqrt((satPos_tTX(i,1) - Xk_m(1))^2 + (satPos_tTX(i,2) - Xk_m(2))^2 + (satPos_tTX(i,3) - Xk_m(3))^2) + Xk_m(4);
end

INOV = pseudorange - h';

Xk_p(:,1) = Xk_m + Kk * (pseudorange - h');


%Update the error covariance matrix
I = eye(5);
Pk_p = (I - Kk * Hk) * Pk_m * (I - Kk * Hk)' + Kk * Rk * Kk';

%Predict the state vector
Xk_pred(:,1) = PHI * Xk_p;

%Predict the error covariance matrix
Pk_pred(:,:,1) = PHI * Pk_p * PHI' + Qk;

%Rename the state estimate variable for easier reading
Xk_est = Xk_p;




end

