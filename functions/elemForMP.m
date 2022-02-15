function [eIN] = elemForMP(coord,etpl,mpC,lp)

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
% [eIN] = ELEMFORMP(coord,etpl,mpC,lp)
%--------------------------------------------------------------------------
% Input(s):
% coord - element coordinates (nen,nD)
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

nD   = size(coord,2);                                                       % number of dimensions
nels = size(etpl,1);                                                        % number of elements
Pmin = mpC-lp;                                                              % particle domain extents (lower)
Pmax = mpC+lp;                                                              % particle domain extents (upper)
a    = true(nels,1);                                                        % initialise logical array
for i=1:nD
  ci = coord(:,i);                                                          % nodal coordinates in current i direction
  c  = ci(etpl);                                                            % reshaped element coordinates in current i direction
  Cmin = min(c,[],2);                                                       % element lower coordinate limit 
  Cmax = max(c,[],2);                                                       % element upper coordainte limit  
  a = a.*((Cmin<Pmax(i)).*(Cmax>Pmin(i)));                                  % element overlap with mp domain
end
eIN = (1:nels).';                                                           % list of all elements
eIN = eIN(a>0);                                                             % remove those elements not in the domain
end                                                                        