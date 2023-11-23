function [mpData] = split_mps(mesh,mpData,directions)
dir_lin = directions(directions ~= 3);
mpData = split_mps_linear(mesh,mpData,dir_lin);
dir_x = zeros(length(mpData),1);
dir_x(1:length(directions)) = directions == 3;
[mpData,new_mps] = split_mps_linear(mesh,mpData,dir_x);
dir_y = zeros(length(mpData),1);
dir_y(new_mps) = 2;
[mpData,new_mps] = split_mps_linear(mesh,mpData,dir_y);
end

