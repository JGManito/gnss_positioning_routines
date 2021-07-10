function [out] = rad2semicircles(in)
%RAD2SEMICIRCLES Converts units from radians into semicircles

out = rad2deg(in) .* ( 1 ./ 180 );

end

