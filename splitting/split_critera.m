function [to_split] = split_critera(mesh,mpData,factor)
%SPLIT_CRITERA Summary of this function goes here
%   Detailed explanation goes here
to_split = zeros(length(mpData),1);
nmp = length(mpData);
bm1=[1 1 1 0 0 0]';
parfor mp=1:nmp 
to_split(mp) = to_split(mp) + (1 * (mpData(mp).lp(1) > (0.5*factor*mesh.h(1))));
to_split(mp) = to_split(mp) + (2 * (mpData(mp).lp(2) > (0.5*factor*mesh.h(2))));
strain = mpData(mp).epsEn;
s=strain-sum(strain(1:3))/3*bm1; j2=(s'*s+s(4:6)'*s(4:6))/2;
if s > 1e-3
    to_split(mp) = 3;
end
%to_split(mp) = 3 * (mpData(mp).lp(1) > (0.5*factor*mesh.h(1))) && (mpData(mp).lp(2) > (0.5*factor*mesh.h(2)));
%to_split(mp) = (mpData(mp).lp(2) > (0.5*factor*mesh.h(2)));
%to_split(mp) = to_split(mp) || (mpData(mp).vp <= 0);
end

