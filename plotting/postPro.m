%Post processing script for the AMPLE code
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% The script produces VTK output files based on the background mesh and
% material point data.  
%
%--------------------------------------------------------------------------
% POSTPRO
%--------------------------------------------------------------------------
% See also:
% MAKEVTK           - VTK file for background mesh
% MAKEVTKMP         - VTK file for MP data
%--------------------------------------------------------------------------

mpDataName = sprintf('output/mpData_%i.vtk',lstp);                          % MP output data file name
sig = reshape([mpData.sig],6,nmp)';                                         % all material point stresses (nmp,6)
mpC = reshape([mpData.mpC],nD,nmp)';                                        % all material point coordinates (nmp,nD)
mpU = [mpData.u]';                                                          % all material point displacements
makeVtkMP(mpC,sig,mpU,mpDataName);                                          % generate material point VTK file
if lstp==1
    makeVtk(mesh.coord,mesh.etpl,'output/mesh.vtk')                         % generate mesh VTK file
end