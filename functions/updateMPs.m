function [mpData] = updateMPs(uvw,mpData)

%Material point update: stress, position and volume
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to update the material point positions and volumes (and domain 
% lengths for GIMPM).  The function also updates the previously converged 
% value of the deformation gradient and the logarithmic elastic strain at 
% each material point based on the converged value and calculates the total 
% displacement of each material point.  
%
% For the generalised interpolation material point method the domain
% lengths are updated according to the stretch tensor following the
% approach of:
% Charlton, T.J., Coombs, W.M. & Augarde, C.E. (2017). iGIMP: An implicit 
% generalised interpolation material point method for large deformations. 
% Computers and Structures 190: 108-125.
%
%--------------------------------------------------------------------------
% [mpData] = UPDATEMPS(uvw,mpData)
%--------------------------------------------------------------------------
% Input(s):
% uvw    - nodal displacements (nodes*nD,1)
% mpData - material point structured array.  The function requires:
%           - mpC : material point coordinates
%           - Svp : basis functions
%           - F   : deformation gradient
%           - lp0 : initial domain lenghts (GIMPM only)
%--------------------------------------------------------------------------
% Ouput(s);
% mpData - material point structured array.  The function modifies:
%           - mpC   : material point coordinates
%           - vp    : material point volume
%           - epsEn : converged elastic strain
%           - Fn    : converged deformation gradient
%           - u     : material point total displacement
%           - lp    : domain lengths (GIMPM only)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------
t = [1 5 9];                                                                % stretch components for domain updating
nmp = length(mpData);                                                       % number of material points
nD  = length(mpData(1).mpC);                                                % number of dimensions
for mp=1:nmp
    nIN = mpData(mp).nIN;                                                   % nodes associated with material point
    nn  = length(nIN);                                                      % number nodes
    N   = mpData(mp).Svp;                                                   % basis functions
    F   = mpData(mp).F;                                                     % deformation gradient                                      
    ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom
    mpU = N*uvw(ed);                                                        % material point displacement
    mpData(mp).mpC   = mpData(mp).mpC + mpU;                                % update material point coordinates
    mpData(mp).vp    = det(F)*mpData(mp).vp0;                               % update material point volumes
    mpData(mp).epsEn = mpData(mp).epsE;                                     % update material point elastic strains
    mpData(mp).Fn    = mpData(mp).F;                                        % update material point deformation gradients
    mpData(mp).u     = mpData(mp).u + mpU.';                                % update material point displacements
    if mpData(mp).mpType == 2                                               % GIMPM only (update domain lengths)        
        [V,D] = eig(F.'*F);                                                 % eigen values and vectors F'F
        U     = V*sqrt(D)*V.';                                              % material stretch matrix        
        mpData(mp).lp = (mpData(mp).lp0).*U(t(1:nD));                       % update domain lengths
    end
end
end