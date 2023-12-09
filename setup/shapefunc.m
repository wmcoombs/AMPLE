function [N,dN] = shapefunc(nen,GpLoc,nD)

%Finite element basis functions
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to provide finite element shapefunctions in 1D, 2D and 3D.  The
% function includes the following elements:
% nen = 8, nD = 3 : tri-linear eight noded hexahedral
% nen = 4, nD = 2 : bi-linear four noded quadrilateral
% nen = 2, nD = 1 : linear two noded line
%
% The function is vectorised so will return the basis functions for any
% number of points.
%
%--------------------------------------------------------------------------
% [N,dN] = SHAPEFUNC(nen,GpLoc,nD)
%--------------------------------------------------------------------------
% Input(s):
% nen    - number of nodes associated with the element
% GpLoc  - local position within the finite element (n,nD)
% nD     - number of dimensions
%--------------------------------------------------------------------------
% Ouput(s);
% N      - shape function matrix (n,nen)
% dN     - derivative of the shape functions (nD,nen)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

n=size(GpLoc,1);                                                            % number of points
N=zeros(n,nen);                                                             % zero shape function matrix
if nD==3                                                                    % 3D
  xsi=GpLoc(:,1); eta=GpLoc(:,2); zet=GpLoc(:,3);                           % local position
  if nen==8                                                                 % 8-noded hexahedral
    N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;
    N(:,2)=(1-xsi).*(1-eta).*(1+zet)/8;
    N(:,3)=(1+xsi).*(1-eta).*(1+zet)/8;
    N(:,4)=(1+xsi).*(1-eta).*(1-zet)/8;
    N(:,5)=(1-xsi).*(1+eta).*(1-zet)/8;
    N(:,6)=(1-xsi).*(1+eta).*(1+zet)/8;
    N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
    N(:,8)=(1+xsi).*(1+eta).*(1-zet)/8;  
  end
elseif nD==2                                                                % 2D
  xsi=GpLoc(:,1); eta=GpLoc(:,2);                                           % local positions
  if nen==4                                                                 % 4-noded quadrilateral
    N(:,1)=0.25*(1-xsi).*(1-eta);
    N(:,2)=0.25*(1-xsi).*(1+eta);
    N(:,3)=0.25*(1+xsi).*(1+eta);
    N(:,4)=0.25*(1+xsi).*(1-eta);
  end
else                                                                        % 1D
  xsi=GpLoc;                                                                % local positions
  if nen==2                                                                 % 2-noded line
    N(:,1)=0.5*(1-xsi);
    N(:,2)=0.5*(1+xsi);  
  end
end


if nargout>1
    dN = zeros(nD,nen);
    if nD==3                                                                    % 3D
        if nen==8                                                                 % 8-noded hexahedral
            dN(1,1)=-(1-eta).*(1-zet)/8;
            dN(1,2)=-(1-eta).*(1+zet)/8;
            dN(1,3)= (1-eta).*(1+zet)/8;
            dN(1,4)= (1-eta).*(1-zet)/8;
            dN(1,5)=-(1+eta).*(1-zet)/8;
            dN(1,6)=-(1+eta).*(1+zet)/8;
            dN(1,7)= (1+eta).*(1+zet)/8;
            dN(1,8)= (1+eta).*(1-zet)/8;
            
            dN(2,1)=-(1-xsi).*(1-zet)/8;
            dN(2,2)=-(1-xsi).*(1+zet)/8;
            dN(2,3)=-(1+xsi).*(1+zet)/8;
            dN(2,4)=-(1+xsi).*(1-zet)/8;
            dN(2,5)= (1-xsi).*(1-zet)/8;
            dN(2,6)= (1-xsi).*(1+zet)/8;
            dN(2,7)= (1+xsi).*(1+zet)/8;
            dN(2,8)= (1+xsi).*(1-zet)/8;
            
            dN(3,1)=-(1-xsi).*(1-eta)/8;
            dN(3,2)= (1-xsi).*(1-eta)/8;
            dN(3,3)= (1+xsi).*(1-eta)/8;
            dN(3,4)=-(1+xsi).*(1-eta)/8;
            dN(3,5)=-(1-xsi).*(1+eta)/8;
            dN(3,6)= (1-xsi).*(1+eta)/8;
            dN(3,7)= (1+xsi).*(1+eta)/8;
            dN(3,8)=-(1+xsi).*(1+eta)/8;
        end
    elseif nD==2                                                                % 2D
        if nen==4                                                                 % 4-noded quadrilateral
            dN(1,1)=-0.25*(1-eta);
            dN(1,2)=-0.25*(1+eta);
            dN(1,3)= 0.25*(1+eta);
            dN(1,4)= 0.25*(1-eta);
            
            dN(2,1)=-0.25*(1-xsi);
            dN(2,2)= 0.25*(1-xsi);
            dN(2,3)= 0.25*(1+xsi);
            dN(2,4)=-0.25*(1+xsi);
        end
    else                                                                        % 1D
        if nen==2                                                                 % 2-noded line
            dN(1)=-0.5;
            dN(2)= 0.5;
        end
    end
end

end