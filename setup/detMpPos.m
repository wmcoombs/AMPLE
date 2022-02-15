function [mpPos] = detMpPos(mp,nD)

%Material point local positions for point generation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to return the local positions of the material points for initial
% set up of a problem.  The material points are evenly spaced through the
% elements that can be regular line, quadrilateral and hexahedral elements
% in 1D, 2D and 3D. 
%
%--------------------------------------------------------------------------
% [mpPos] = DETMPPOS(mp,nD)
%--------------------------------------------------------------------------
% Input(s):
% mp    - number of material points in each direction
% nD    - number of dimensions
%--------------------------------------------------------------------------
% Ouput(s);
% mpPos - local material point positions (nmp,nD)
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

nmp=mp^nD;                                                                  % number of material points per element
mpPos=zeros(nmp,nD);                                                        % zero material point positions
a=2/mp;                                                                     % local length associated with the material point
b=(a/2:a:2)-1;                                                              % local positions of MP in 1D
if nD==1                                                                    % 1D
  mpPos=b.';                                                            
elseif nD==2                                                                % 2D
  for i=1:mp
    for j=1:mp
      mpPos((i-1)*mp+j,1)=b(i);
      mpPos((i-1)*mp+j,2)=b(j);
    end
  end
else                                                                        % 3D
  for i=1:mp
    for j=1:mp
      for k=1:mp 
        mpPos((i-1)*mp^2+(j-1)*mp+k,1)=b(i);
        mpPos((i-1)*mp^2+(j-1)*mp+k,2)=b(j);
        mpPos((i-1)*mp^2+(j-1)*mp+k,3)=b(k);
      end 
    end
  end  
end
end