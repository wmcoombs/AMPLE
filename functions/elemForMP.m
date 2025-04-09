function [eIN] = elemForMP(mesh,mpC,lp)

%Find elements associated with the material point
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   06/05/2015
% Description:
% Function to determine the elements that are associated with a material
% point assuming that the material point's domain is symmetric about the
% particle position.
%
%--------------------------------------------------------------------------
% [eIN] = ELEMFORMP(mesh,mpC,lp)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh structured array. Function requires
%           - etpl  : element topology (nels,nen) 
%           - eMin  : element lower coordinate limit (nels,nD)
%           - eMax  : element upper coordinate limit (nels,nD)
% etpl  - element topology (nels,nen)
% mpC   - material point coordinates (1,nD)
% lp    - domain half width
%--------------------------------------------------------------------------
% Ouput(s);
% eIN   - vector containing the elements associated with the mp
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

Cmin = mesh.eMin;                                                           % element lower coordinate limit 
Cmax = mesh.eMax;                                                           % element upper coordinate limit 
etpl = mesh.etpl;                                                           % element topology
nD   = size(Cmin,2);                                                        % number of dimensions
nels = size(etpl,1);                                                        % number of elements
d    = sqrt(eps)*lp;                                                        % tolerance on the search (avoids zero overlap detection with GIMPM)
Pmin = mpC-lp+d;                                                            % particle domain extents (lower)
Pmax = mpC+lp-d;                                                            % particle domain extents (upper)
a    = true(nels,1);                                                        % initialise logical array
for i=1:nD 
    a = a.*((Cmin(:,i)<Pmax(i)).*(Cmax(:,i)>Pmin(i)));                      % element overlap with mp domain
end
eIN = find(a);                                                              % elements overlaps by the domain
end                                                                        