function [Pk_prev,Xk_est,Xk_pred] = positionDataKF_DD(flag,sigmaR,P_initial,measurement,H,x0,Pk_prev,Xk_prev)
%POSITIONDATAEKF This function implements a Kalman Filter for use
%with receiver position data.
%   Detailed explanation goes here

if strcmp('initialize',flag) == 1
    n = size(measurement,1);
    Xk_prev = x0;                   %Initial position estimate
    Pk_prev = P_initial .* eye(3);    %Initial error covariance matrix
    PHI = eye(3);                   %State transition matrix
    Q = zeros(3,3);                 %Noise covariance matrix of the continuous time model
    R = sigmaR^2 * eye(n);          %Observation noise covariance matrix
    zk = measurement;               %Measurement vector
    
    
else
    n = size(measurement,1);
    PHI = eye(3);               %State transition matrix
    Q = zeros(3,3);             %Noise covariance matrix of the continuous time model
    R = sigmaR^2 * eye(n);      %Observation noise covariance matrix
    zk = measurement;           %Measurement vector
end
    
    %Convert line vectors to column vectors
    Xk_prev = Xk_prev';
    %zk = zk';
    
    %Prediction step
    Xk_m = PHI * Xk_prev;       %Predict the state vector
    Pk_m = PHI * Pk_prev * PHI' + Q; %Predict the error covariance matrix
    
    %Compute the Kalman Gain
    Kk = Pk_m * H' * inv(H * Pk_m * H' + R);
    
    %Compute the state estimate
    Xk_p = Xk_m + Kk * (zk - H*Xk_m);
    
    %Compute the error covariance
    Pk_p = Pk_m - Kk * H * Pk_m;
    
    %Rename the output variables
    Xk_est = Xk_p';
    Xk_pred = Xk_m';
    Pk_prev = Pk_p;
    
    
    
    
end

