function [mpData] = updateMPs(uvw,mpData,mesh)

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
    if mpData(mp).vp <= 0
        throw(MException("AMPLE:negative_volume","MP volume is negative\n"));
    end
    if mpData(mp).mpType == 2                                               % GIMPM only (update domain lengths)   
        corner_tracking = false;
        if corner_tracking == false
        [V,D] = eig(F.'*F);                                                 % eigen values and vectors F'F
        U     = V*sqrt(D)*V.';                                              % material stretch matrix        
        mpData(mp).lp = (mpData(mp).lp0).*U(t(1:nD));                       % update domain lengths
        else
        for c = 1:4
            N = mpData(mp).CSvp(c).Svp;
            nIN = mpData(mp).CNodes(c).N;
            nn  = length(nIN);
            ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom
            mpU = N*uvw(ed);
            mpData(mp).C(c,:) = mpData(mp).C(c,:) + mpU;
        end
        lp = [0.5*(max(mpData(mp).C(:,1)) - min(mpData(mp).C(:,1)));
              0.5*(max(mpData(mp).C(:,2)) - min(mpData(mp).C(:,2)))];
        nom = det(F(1:nD,1:nD))*mpData(mp).lp0(1)*mpData(mp).lp0(2);
        den = lp(1)*lp(2);
        mpData(mp).lp(1) = lp(1)*(nom/den).^(1/2);
        mpData(mp).lp(2) = lp(2)*(nom/den).^(1/2);
        x = mpData(mp).lp(1);
        y = mpData(mp).lp(2);

        % 
        if ((-x + mpData(mp).mpC(1))<=0)
            x = mpData(mp).mpC(1)-10^-7;
        end
        max_X = mesh.size(1);
        if ((x + mpData(mp).mpC(1))>max_X)
            x = (max_X - mpData(mp).mpC(1))-10^-7;
        end
        mpData(mp).lp(1) = x;


        if ((-y + mpData(mp).mpC(2))<=0)
            y = mpData(mp).mpC(2)-10^-7;
        end
        max_Y = mesh.size(2);
        if ((y + mpData(mp).mpC(2))>max_Y)
            y = (max_Y - mpData(mp).mpC(2))-10^-7;
        end
        mpData(mp).lp(2) = y;
        mpData(mp).C = mpCorners(x,y) + repmat(mpData(mp).mpC,4,1);
        end
    end
end
end