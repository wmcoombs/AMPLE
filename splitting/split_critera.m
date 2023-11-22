function [to_split] = split_critera(mesh,mpData,factor)
%SPLIT_CRITERA Summary of this function goes here
%   Detailed explanation goes here
to_split = zeros(length(mpData),1);
nmp = length(mpData);
parfor mp=1:nmp 
%to_split(mp) = (mpData(mp).lp(1) > (0.5*factor*mesh.h(1))) || (mpData(mp).lp(2) > (0.5*factor*mesh.h(2)));
to_split(mp) = (mpData(mp).lp(2) > (0.5*factor*mesh.h(2)));
%to_split(mp) = to_split(mp) || (mpData(mp).vp <= 0);
end

