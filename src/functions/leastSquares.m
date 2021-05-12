function [x0,dx,GDOP,PDOP,TDOP,H] = leastSquares(observation,satPos_tTX,initialEstimate)
%LEASTSQUARES Summary of this function goes here
%   Detailed explanation goes here

%Define constants
c = 2.99792458 * 10^8; %WGS-84 definition of speed of light

%Define the required precision
prec = 10^-3;

%Split the input arrays for easier reading
pseudorange = observation(:,5);

%Define initial position estimate
x0 = [initialEstimate(1),initialEstimate(2),initialEstimate(3),0]; %Center of the Earth in an ECEF frame


%Keep computing the Least Squares solution until the error is below the
%required precision

posError = [1,1,1];
while norm(posError) >= prec
    
    %Compute the geometric range between the initial receiver position estimate
    %and the satellite positions
    for i = 1:size(observation,1)
        rho0(i) = norm(satPos_tTX(i,:) - x0(1:3));
    end
    
    %Compute the drho vector
    drho = []; %Clear previous result
    for i = 1:size(observation,1)
        drho(i) = pseudorange(i) - (rho0(i) + x0(4));
    end
    drho = transpose(drho);
    
    %Compute the elements of matrix H
    for i = 1:size(observation,1)
        H(i,1) = (x0(1)-satPos_tTX(i,1))/(rho0(i));
        H(i,2) = (x0(2)-satPos_tTX(i,2))/(rho0(i));
        H(i,3) = (x0(3)-satPos_tTX(i,3))/(rho0(i));
        H(i,4) = 1;
    end
    
    %Compute the elements of the dx vector using the least squares method
    dx = inv(H'*H) * H' *drho;
    %Update the position solution
    x0 = x0 + dx';
    posError = dx(1:3);
    
    
end
M = inv(H'*H);
GDOP = sqrt(trace(M));
PDOP = sqrt(trace(M(1:3,1:3)));
TDOP = sqrt(M(4,4));

end

