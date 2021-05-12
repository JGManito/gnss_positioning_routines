function [duration,meanX,meanY,meanZ,meanV,nObs,valid,svinActive] = parseSVIN(inputData)
%PARSESVIN Summary of this function goes here
%   Detailed explanation goes here

count = 0;
i = 1;
ii = 1;
for i = 1:size(inputData,1)
    if inputData(i) == hex2dec('B5')
        if inputData(i+1) == hex2dec('62') %Start of a message
            if inputData(i+2) == hex2dec('0D')
                if inputData(i+3) == hex2dec('04') %TIM-SVIN message
                    if double(inputData(i+31)) == 1  && double(inputData(i+30)) ~= 1 %If SVIN active and no valid flag
                        count = count +1;
                        duration(ii) = readDataFormat([inputData(i+6),inputData(i+7),inputData(i+8),inputData(i+9)],'U4');
                        
                        meanX(ii) = readDataFormat([inputData(i+10),inputData(i+11),inputData(i+12),inputData(i+13)],'I4') * 0.01; %Data in centimeters
                        meanY(ii) = readDataFormat([inputData(i+14),inputData(i+15),inputData(i+16),inputData(i+17)],'I4') * 0.01; %Data in centimeters
                        meanZ(ii) = readDataFormat([inputData(i+18),inputData(i+19),inputData(i+20),inputData(i+21)],'I4') * 0.01; %Data in centimeters
                        
                        meanV(ii) = readDataFormat([inputData(i+22),inputData(i+23),inputData(i+24),inputData(i+25)],'U4') * 0.000001; %Data in mm^2
                        
                        nObs(ii) = readDataFormat([inputData(i+26),inputData(i+27),inputData(i+28),inputData(i+29)],'U4');
                        
                        valid(ii) = readDataFormat(inputData(i+30),'U1');
                        
                        svinActive(ii) = readDataFormat(inputData(i+31),'U1');
                        
                        svinActive(ii) = double(inputData(i+31));
                        
                        ii = ii+1;
                    end
                    
                    if valid(ii) == 1
                        break;
                    end
                end
            end
        end
    end
    
    
    
end

function out = readDataFormat (in,dataType)

switch dataType
    case 'U4'
        out = double(typecast([uint8(in(1)),uint8(in(2)),uint8(in(3)),uint8(in(4))],'uint32'));
        
    case 'I4'
        out = double(typecast([uint8(in(1)),uint8(in(2)),uint8(in(3)),uint8(in(4))],'int32'));
        
    case 'U1'
        out = double(in);
        
        
end