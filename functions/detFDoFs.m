function [fd] = detFDoFs(mesh)
%Determine the free degrees of freedom on the background mesh
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   17/12/2018
% Description:
% Function to determine the free degrees of freedom of the background mesh
% based on the elements that contain material points and the displacement
% boundary conditions. 
%
%--------------------------------------------------------------------------
% [fd] = DETFDOFS(etpl,eInA,bc,nD,nDoF)
%--------------------------------------------------------------------------
% Input(s):
% mesh  - mesh structured array. Function requires: 
%           - etpl : element topology (nels,nen) 
%           - eInA : elements "active" in the analysis
%           - bc   : boundary conditions (*,2)
%--------------------------------------------------------------------------
% Ouput(s);
% fd    - free degrees of freedom on the background mesh (*,1)
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

[nodes,nD] = size(mesh.coord);                                              % no. nodes and dimensions
nDoF   = nodes*nD;                                                          % no. degrees of freedom
incN   = unique(mesh.etpl(mesh.eInA>0,:));                                  % unique active node list
iN     = size(incN,1);                                                      % number of nodes in the list
incDoF = reshape(ones(nD,1)*incN'*nD-(nD-1:-1:0).'*ones(1,iN),1,iN*nD);     % active degrees of freedom
fd     = (1:nDoF);                                                          % all degrees of freedom
fd(mesh.bc(:,1)) = 0;                                                       % zero fixed displacement BCs
fd     = fd(incDoF);                                                        % only include active DoF 
fd     = fd(fd>0);                                                          % remove fixed displacement BCs
end