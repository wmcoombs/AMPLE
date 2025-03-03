%AMPLE GS: A Material Point Learning Environment with Ghost Stabilisation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   03/03/2025
% Description:
% Large deformation elasto-plastic (EP) material point method (MPM) code
% based on an updated Lagrangian (UL) descripition of motion with a 
% quadrilateral background mesh. The implementation includes ghost
% stabilisation based on https://doi.org/10.1002/nme.7332
%
%--------------------------------------------------------------------------
% See also:
% SETUPGRID             - analysis specific information
% ELEMMPINFO            - material point-element information
% GHOSTPENALTY          - Ghost stabilisation 
% DETEXTFORCE           - external forces
% DETFDOFS              - mesh unknown degrees of freedom
% LINSOLVE              - linear solver
% DETMPS                - material point stiffness and internal force
% UPDATEMPS             - update material points
% POSTPRO               - post processing function including vtk output
%--------------------------------------------------------------------------
clear;
addpath('constitutive','functions','plotting','setup');        
[lstps,g,mpData,mesh] = setupGrid;                                          % setup information
NRitMax = 10; tol = 1e-6;                                                   % Newton Raphson parameters
[nodes,nD] = size(mesh.coord);                                              % number of nodes and dimensions
[nels,nen] = size(mesh.etpl);                                               % number of elements and nodes/element
nDoF = nodes*nD;                                                            % total number of degrees of freedom
nmp  = length(mpData);                                                      % number of material points
lstp = 0;                                                                   % zero loadstep counter (for plotting function)
uvw  = zeros(nDoF,1);                                                       % zeros displacements (for plotting function)
run postPro;                                                                % plotting initial state & mesh
for lstp=1:lstps                                                            % loadstep loop
  fprintf(1,'\n%s %4i %s %4i\n','loadstep ',lstp,' of ',lstps);             % text output to screen (loadstep)
  [mesh,mpData] = elemMPinfo(mesh,mpData);                                  % material point - element information
  Ks   = mesh.ks*ghostPenalty(mesh);                                        % Ghost stabilisation matrix
  fext = detExtForce(nodes,nD,g,mpData);                                    % external force calculation (total)
  fext = fext*lstp/lstps;                                                   % current external force value
  oobf = fext;                                                              % initial out-of-balance force
  fErr = 1;                                                                 % initial error
  frct = zeros(nDoF,1);                                                     % zero the reaction forces
  uvw  = zeros(nDoF,1);                                                     % zero the displacements
  fd   = detFDoFs(mesh);                                                    % free degrees of freedom
  NRit = 0;                                                                 % zero the iteration counter
  Kt   = sparse(nDoF,nDoF);                                                 % zero global stiffness matrix
  while (fErr > tol) && (NRit < NRitMax) || (NRit < 2)                      % global equilibrium loop
    [duvw,drct] = linSolve(mesh.bc,Kt+Ks,oobf,NRit,fd);                     % linear solver
    uvw  = uvw+duvw;                                                        % update displacements
    frct = frct+drct;                                                       % update reaction forces
    fgst = Ks*uvw;                                                          % ghost stabilisation force
    [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
    oobf = (fext-fint+frct-fgst);                                           % out-of-balance force
    fErr = norm(oobf)/norm(fext+frct+eps);                                  % normalised oobf error
    NRit = NRit+1;                                                          % increment the NR counter
    fprintf(1,'%s %2i %s %8.3e\n','  iteration ',NRit,' NR error ',fErr);   % text output to screen (NR error)
  end                           
  mpData = updateMPs(uvw,mpData);                                           % update material points
  run postPro;                                                              % Plotting and post processing   
end


