function [float_value] = exp2float(exp_value)
%EXP2FLOAT This function reads a string in the Dn.m format and outputs it
%as a double-precision floating point number (according to IEEE® Standard 754)
%   Detailed explanation goes here

[exp_split,~] = split(exp_value,["D","E","e"]);

frac = exp_split{1};
exp = exp_split{2};

frac = str2num(frac);
exp = str2num(exp);

float_value = frac * 10^(exp);


end

