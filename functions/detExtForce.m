function [fext] = detExtForce(nodes,nD,g,mpData)

%Global external force determination  
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to determine the external forces at nodes based on body forces
% and point forces at material points.
%
%--------------------------------------------------------------------------
% [fbdy,mpData] = DETEXTFORCE(coord,etpl,g,eIN,mpData)
%--------------------------------------------------------------------------
% Input(s):
% nodes  - number of nodes (total in mesh)
% nD     - number of dimensions
% g      - gravity
% mpData - material point structured array. Function requires:
%           mpM   : material point mass
%           nIN   : nodes linked to the material point
%           Svp   : basis functions for the material point
%           fp    : point forces at material points
%--------------------------------------------------------------------------
% Ouput(s);
% fext   - external force vector (nodes*nD,1)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

nmp  = size(mpData,2);                                                      % number of material points & dimensions 
fext = zeros(nodes*nD,1);                                                   % zero the external force vector
grav = zeros(nD,1); grav(nD) = -g;                                          % gavity vector
for mp = 1:nmp
   nIN = mpData(mp).nIN;                                                    % nodes associated with MP
   nn  = length(nIN);                                                       % number of nodes influencing the MP
   ed  = reshape(ones(nD,1)*(nIN-1)*nD+(1:nD).'*ones(1,nn),1,nn*nD);        % node degrees of freedom
   Svp = mpData(mp).Svp;                                                    % basis functions
   fext(ed) = fext(ed)+mpData(mp).mpM*reshape(grav*Svp,nn*nD,1)...          % global body force contribution
            + reshape(mpData(mp).fp*Svp,nn*nD,1);                           % material point force contribution
end
