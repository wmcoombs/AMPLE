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

[~,nD]=size(coord); [nels,nen]=size(etpl);                                  % basic size information
Pmin=mpC-lp; Pmax=mpC+lp;                                                   % particle domain extents
a=true(nels,1);
for i=1:nD
  c=reshape(coord(reshape(etpl',nels*nen,1),i),nen,nels).';                 % reshaped element coordinates
  Cmin=min(c.').';                                                          % element lower coordinate limit
  Cmax=max(c.').';                                                          % element upper coordainte limit  
  a=a.*((Cmin<Pmax(i)).*(Cmax>Pmin(i)));                                    % element overlap with mp domain
end
eIN=(1:nels).';                                                             % list of all elements
eIN=eIN(a>0);                                                               % remove those elements not in the domain