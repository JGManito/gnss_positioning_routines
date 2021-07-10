function [floatOut,remainderStr] = splitDataFormat(inputStr,format)
%SPLITDATAFORMAT This function receives a string of data to compare against
%a known data format and returns the value and the rest of the string to
%allow chaining of this function multiple times
%   Possible formats:
%       In - Integer with n digits
%       Fn.m - Float with n total characters (including signs and decimal
%       point) with m decimal digits

%Choose method based on the type of variable

switch format(1)
    case 'I' %Integer data format
        n = str2num(format(2:end)); %Get the number of digits
        
        strAcc = []; %Create a string to work as an accumulator for the parsing of the number
        for i = 1:n
            if inputStr(i) == ' ' %If empty space
                strAcc = strcat(strAcc,'0');
            else
                strAcc = strcat(strAcc,inputStr(i));
            end
        end
        floatOut = str2num(strAcc);
        remainderStr = inputStr(n+1:end);
        
        %-----------------------------------------------%
    case 'F' %Float data format
        digits = strsplit(format(2:end),'.'); %Get the number of total and decimal digits
        n = str2num(digits{1});
        m = str2num(digits{2}); %Possibly not necessary due to using a string as accumulator
        
        strAcc = []; %Create a string to work as an accumulator for the parsing of the number
        for i = 1:n
            if inputStr(i) == ' ' %If empty space
                strAcc = strcat(strAcc,'0');
            else
                strAcc = strcat(strAcc,inputStr(i));
            end
        end
        floatOut = str2num(strAcc);
        remainderStr = inputStr(n+1:end);
end

end

