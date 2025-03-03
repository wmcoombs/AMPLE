function [GpLoc] = findLocalCoord(xGp,coord)

%Determine local coordinates
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/06/2021
% Description:
% Function to determine the local [position of a point within a finite
% element based on knowing the global position wihtin the element.  
% 
%--------------------------------------------------------------------------
% [GpLoc] = FINDLOCALCOORD(xGp,coord)
%--------------------------------------------------------------------------
% Input(s):
% xGp    - global position of the point in the element (1,D)
% coord  - global coordinates of the element ordered according to the local
%          convention (nen,nD)
%--------------------------------------------------------------------------
% Ouput(s);
% GpLoc  - local position in the element (1,nD)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

tol=1e-9;                                                                   % convergence tolerance
[nen,nD]=size(coord);                                                       % no. element nodes and dimensions
GpLoc=zeros(1,nD);                                                          % zero position in element
[N,dNr] = shapefunc(nen,GpLoc,nD);                                          % shape functions & local derivatives
while norm(xGp-N*coord)>tol                                                 % Newton loop
  GpLoc = GpLoc+(xGp-N*coord)/(dNr*coord);                                  % local position in element
  [N,dNr] = shapefunc(nen,GpLoc,nD);                                        % shape functions & local derivatives
end

end