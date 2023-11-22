function [Svp,dSvp] = MPMbasis(mesh,mpData,node,cmat)

%Basis functions for the material point method
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to determine the multi-dimensional MPM shape functions from the
% one dimensional MPM functions.  The function includes both the standard
% and generalised interpolation material point methods.
%
%--------------------------------------------------------------------------
% [Svp,dSvp] = MPMBASIS(coord,mpC,L)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh data structured array. Function requires:
%           - coord  : nodal coordinates
%           - h      : grid spacing
%
% mpData - material point structured array.  Function requires:
%           - mpC    : material point coordinates (single point)
%           - lp     : particle domain lengths
%           - mpType : material point type (1 or 2)
%
% node   - background mesh node number
%--------------------------------------------------------------------------
% Ouput(s);
% Svp   - particle characteristic function
% dSvp  - gradient of the characterstic function
%--------------------------------------------------------------------------
% See also:
% SVPMPM    - MPM basis functions in 1D (mpType = 1
% SVPGIMP   - GIMPM basis functions in 1D (mpType = 2)
%--------------------------------------------------------------------------

coord  = mesh.coord(node,:);                                                % node coordinates
mpC    = mpData.mpC;                                                        % material point coordinates
lp     = mpData.lp;                                                         % material point domain length
mpType = mpData.mpType;                                                     % material point type (MPM or GIMPM)
nD     = size(mpC,1)*size(mpC,2);                                           % number of dimensions
h = mesh.h;
if nargin == 4
    mpC = cmat;
    lp = zeros(size(lp));
    mpType = 1;
end

S=zeros(nD,1); dS=S; dSvp=S;                                                % zero vectors used in calcs
for i=1:nD
    if mpType == 1
        h_XY   = h(i);                                % grid spacing
        [S(i),dS(i)] = SvpMPM(mpC(i),coord(i),h_XY);                        % 1D MPM functions
    elseif mpType == 2
        h_XY   = h(i);
        [S(i),dS(i)] = SvpGIMP(mpC(i),coord(i),h_XY,lp(i));                 % 1D GIMPM functions
    end
end
if nD == 1
    indx = [];                                                              % index for basis derivatives (1D)
elseif nD == 2
    indx = [2; 1];                                                          % index for basis derivatives (2D)
elseif nD == 3
    indx = [2 3; 1 3; 1 2];                                                 % index for basis derivatives (3D)
end
Svp=prod(S);                                                                % basis function
for i=1:nD
    dSvp(i)=dS(i)*prod(S(indx(i,:)));                                       % gradient of the basis function
end
end