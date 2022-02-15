function [mesh,mpData] = elemMPinfo(mesh,mpData)

%Determine the basis functions for material points 
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to determine the basis functions and spatial derivatives of each
% material point.  The function works for regular background meshes with
% both the standard and generalised interpolation material point methods.
% The function also determines, and stores, the elements associated with
% the material point and a unique list of nodes that the material point
% influences.  The number of stiffness matrix entries for each material
% point is determined and stored. 
%
%--------------------------------------------------------------------------
% [fbdy,mpData] = ELEMMPINFO(mesh,mpData)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh structured array. Function requires: 
%           - coord : coordinates of the grid nodes (nodes,nD)
%           - etpl  : element topology (nels,nen) 
%           - h     : background mesh size (nD,1)
% mpData - material point structured array.  Function requires:
%           - mpC   : material point coordinates
%--------------------------------------------------------------------------
% Ouput(s);
% mesh   - mesh structured array. Function modifies:
%           - eInA  : elements in the analysis 
% mpData - material point structured array. Function modifies:
%           - nIN   : nodes linked to the material point
%           - eIN   : element associated with the material point
%           - Svp   : basis functions for the material point
%           - dSvp  : basis function derivatives (at start of lstp)
%           - nSMe  : number stiffness matrix entries for the MP
%--------------------------------------------------------------------------
% See also:
% ELEMFORMP         - find elements for material point
% NODESFORMP        - nodes associated with a material point 
% MPMBASIS          - MPM basis functions
%--------------------------------------------------------------------------

nmp      = size(mpData,2);                                                  % number of material points & dimensions
[~,nD]   = size(mesh.coord);                                                % number of nodes (total in mesh)  
[nels,~] = size(mesh.etpl);                                                 % numerb of elements in mesh
mpC  = reshape([mpData.mpC],nD,nmp).';                                      % all material point coordinates (nmp,nD)
lp   = reshape([mpData.lp] ,nD,nmp).';                                      % all domain lengths
eInA = zeros(nels,1);                                                       % zero elements taking part in the analysis
for mp = 1:nmp
    eIN  = elemForMP(mesh.coord,mesh.etpl,mpC(mp,:),lp(mp,:));              % elements connected to the material point
    nIN  = nodesForMP(mesh.etpl,eIN).';                                     % unique list of nodes associated with elements
    nn   = length(nIN);                                                     % number of nodes influencing the MP
    Svp  = zeros(1,nn);                                                     % zero basis functions
    dSvp = zeros(nD,nn);                                                    % zero basis function derivatives
    for i = 1:nn
        node = nIN(i);                                                      % current node
        [S,dS] = MPMbasis(mesh,mpData(mp),node);                            % basis function and spatial derivatives
        Svp(i) = Svp(i) + S;                                                % basis functions for all nodes
        dSvp(:,i) = dSvp(:,i) + dS;                                         % basis function derivatives
    end  
    mpData(mp).nIN  = nIN;                                                  % nodes associated with material point
    mpData(mp).eIN  = eIN;                                                  % elements associated with material point 
    mpData(mp).Svp  = Svp;                                                  % basis functions
    mpData(mp).dSvp = dSvp;                                                 % basis function derivatives
    mpData(mp).nSMe = (nn*nD)^2;                                            % number stiffness matrix components
    eInA(eIN) = 1;                                                          % identify elements in the analysis
end
mesh.eInA = eInA;                                                           % store eInA to mesh structured array
end