function [varargout] = RAIM(flag,varargin)
%RAIM Summary of this function goes here
%   Detailed explanation goes here

switch flag
    case 'initialize'
        %Only input is the required probability of false alarm, Pfa
        if nargin == 2
            Pfa = varargin{1};
        else
            disp("Too many input variables for RAIM initialization")
            return
        end
        
        %Initialize RAIM
        %Compute the threshold for nsats=5:16 (minimum and maximum number
        %of satellites in view)
        for n=5:16
            %1. Get parameter a
            a = (n - 4)/2;
            
            %2. Compute the gamma function
            gamma_a = gamma(a);
            
            %3. Compute the threshold using the lower incomplete gamma function
            func = @(s) (1./(2^a * gamma_a)).*(exp(-s/2) .* s.^(a-1));
            integralfunc = @(lambdaVar) (integral(func,0,lambdaVar.^2) - 1 + Pfa);
            lambda(n) = fzero(integralfunc,5);
        end
        varargout{1} = lambda;
        
    case {'run','exclude'}
        %Split the input variables
        if nargin == 6
            observation = varargin{1};
            navMessage = varargin{2};
            satPos_tTX = varargin{3};
            threshold = varargin{4};
            initialEstimate = varargin{5};
            
        elseif nargin == 7
            observation = varargin{1};
            navMessage = varargin{2};
            satPos_tTX = varargin{3};
            threshold = varargin{4};
            initialEstimate = varargin{5};
            H = varargin{6};
        else
            disp("Wrong number of input parameters for RAIM algorithm")
            return
        end
        
        %Expand input variables for easier parsing
        pseudoranges = observation(:,5);
        nSats = size(observation,1);
        
        %Compute the least squares solution
        initialEstimate = [0,0,0]; %Center of the Earth in ECEF coordinates.
        if exist('H')
            %Compute just the position solution
            [recPos,drecPos,~,~,~,~] = leastSquares(observation,satPos_tTX,initialEstimate);
        else
            %Least-Squares H matrix not present. Compute it and the
            %position solution
            [recPos,drecPos,~,~,~,H] = leastSquares(observation,satPos_tTX,initialEstimate);
        end
        
        %Get the geometric range between receiver and satellite
        for i = 1:size(observation,1)
            r(i) = norm(satPos_tTX(i,:) - recPos(1:3));
        end
        
        %Get the delta pseudorange
        for i = 1:size(observation,1)
            dpseudorange(i) = pseudoranges(i) - (r(i) + recPos(4));
        end
        dpseudorange = transpose(dpseudorange);
        
        %Get the estimated pseudoranges
        dpseudorangeEst = H*drecPos;
        
        %Compute the pseudorange residue
        w = dpseudorange - dpseudorangeEst;
        
        %Compute the sum of squared errors
        SSE = transpose(w)*w;
        
        %Compute the test statistic
        t = sqrt(SSE/(nSats-4));
        
        if t >= threshold(nSats)
            if strcmp(flag,'exclude')
                %In exclude mode just run FD, we're looking for a
                %combination that doesn't trigger FD.
                outputLevel = 1;
                outputObs = observation;
                outputNavMessage = navMessage;
                outputSatPos = satPos_tTX;
                
                %Apply FD or FDE depending on the number of satellites
            elseif nSats == 5
                outputLevel = 1;
                outputObs = observation;
                outputNavMessage = navMessage;
                outputSatPos = satPos_tTX;
                fprintf("RAIM FD @t=%d\n",observation(1,3));
                
            else
                %To-Do: Find the combination of satellites that minimizes
                %the error!
                outputLevel = 2;
                fprintf("RAIM FDE @t=%d\n",observation(1,3));
                
                %Create an array that lists all the possible combinations
                %of the indexes of the observation array
                combinations = nchoosek(1:nSats,nSats-1);
                
                for i=1:size(combinations,1)
                    %Create an array of observations and navMessages
                    %with only N-1 satellites
                    for j=1:size(combinations,2)
                        index = combinations(i,j);
                        observationTemp(j,:) = observation(index,:);
                        navMessageTemp(j,:) = navMessage(index,:);
                        satPositionTemp(j,:) = satPos_tTX(index,:);
                    end
                    
                    %Run the RAIM algorithm on the subset
                    [FD_outputLevel,~,~,~] = RAIM('exclude',observationTemp,navMessageTemp,satPositionTemp,threshold,initialEstimate);
                    
                    %Remove the satellite that degrades the position
                    %solution
                    if FD_outputLevel == 0
                        
                        %Find the excluded SVN
                        set = 1:nSats;
                        missingIndex = find(ismember(set,combinations(i,:)) == 0,1,'first');
                        excludedSVN = observation(missingIndex,4);
                        fprintf("RAIM FDE @t=%d: Removed satellite %.0f from the constellation\n",observation(1,3),excludedSVN);
                        outputLevel = 2;
                        outputObs = observationTemp;
                        outputNavMessage = navMessageTemp;
                        outputSatPos = satPositionTemp;
                        break;
                    end
                    
                    %If the above condition hasn't been reached, there's
                    %more than 1 problem satellite.
                    if FD_outputLevel == 1 && i == size(combinations,1) 
                        outputLevel = -1;
                        outputObs = observation;
                        outputNavMessage = navMessage;
                        outputSatPos = satPos_tTX;
                    end
                    
                end
            end
        else
            outputLevel = 0;
            outputObs = observation;
            outputNavMessage = navMessage;
            outputSatPos = satPos_tTX;
        end
        
        varargout{1} = outputLevel;
        varargout{2} = outputObs;
        varargout{3} = outputNavMessage;
        varargout{4} = outputSatPos;
        varargout{5} = t;
        
        
        
        
end







end

