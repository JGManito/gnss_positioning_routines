function [cellOut] = splitNum(inputStr)
%SPLITNUM Splits a string with several numerical values joined together without
%whitespaces and splits it
%   This function works similarly to the split function from matlab,
%   however with extra functionality to split numerical values that are joined
%   together. The input is a whole string to parse, and outputs the split
%   values exactly like they are written in the string (no conversion
%   happens)


%Do the initial split with whitespaces as separators
splitStr = split(inputStr);

%Find decimal points and split the string using them as reference (all
%numbers are in scientific notation)
for i=1:size(splitStr,1) %Go through each cell
    aux = splitStr{i}; %Save cell content as string to allow for easier parsing
    
    k=1;
    for j=1:size(aux,2)
        if aux(j) == '.'
            if (j>2 && aux(j-2) == '-') %Check for minus sign
                split_position(i,k) = j - 2; %Create an array of the string positions where to split
                k = k + 1;
            else
                split_position(i,k) = j - 1;
                k = k + 1;
            end
        end
    end
end

k = 1; %Counter for the index of the array that contains the parsed strings

for i=1:size(split_position,1)
    if split_position(i,1) ~= 0 %It's only 0 if there's nothing to split on that cell
        aux = splitStr{i}; %Save cell content as string to allow for easier parsing
        
        strAux = aux(1:split_position(i,1)-1);
        for j=1:size(split_position,2) %Iterate for all the split positions
            if split_position(i,j) ~= 0
                if j<size(split_position,2) && split_position(i,j+1)~=0 %Split before the end of the string
                    strAux = strcat(strAux," ",aux(split_position(i,j):split_position(i,j+1)-1));
                else %Split to the end of the string
                    if isempty(strAux) %In case a number was split correctly initially but still got recognized because it had a decimal point
                        strAux = aux(split_position(i,j):end);
                    else
                        strAux = strcat(strAux," ",aux(split_position(i,j):end));
                    end
                end
            end
        end
        splitStr (i,:) = cellstr(strAux);
        %k = k+1;
    end
end

%Rejoin the cells and split the output again (to avoid having to work
%around the changing of the index numbers in the cell array)
aux = join(splitStr);
cellOut = split(aux);

%Remove empty cells that might have been left as parsing artifacts
cellOut = cellOut(~cellfun('isempty',cellOut));

end

