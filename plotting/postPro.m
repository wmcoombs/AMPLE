%Post processing script for the AMPLE code
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% The script produces VTK output files based on the background mesh and
% material point data.  
% Background mesh is plotted for all loadsteps.  
%
%--------------------------------------------------------------------------
% POSTPRO
%--------------------------------------------------------------------------
% See also:
% MAKEVTK           - VTK file for background mesh
% MAKEVTKMP         - VTK file for MP data
%--------------------------------------------------------------------------

sig = reshape([mpData.sig],6,nmp)';                                         % all material point stresses (nmp,6)
mpC = reshape([mpData.mpC],nD,nmp)';                                        % all material point coordinates (nmp,nD)
mpU = [mpData.u]';                                                          % all material point displacements
lps = [mpData.lp];
lp = [lps(1:2:end);lps(2:2:end)]';
plastic = reshape([mpData.epsPlastic],6,nmp)';
plastic = plastic(:,1);
mpDataName = sprintf('output/mpData_%i.vtk',lstp);                          % MP output data file name
makeVtkMP(mpC,sig,mpU,lp,plastic,mpDataName);                                          % generate material point VTK file

meshName = sprintf('output/mesh_%i.vtk',lstp);                              % MP output data file name
makeVtk(mesh.coord,mesh.etpl,uvw,meshName);                                 % generate mesh VTK file