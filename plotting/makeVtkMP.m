function makeVtkMP(mpC,sig,uvw,lp,plastic,mpFileName)

%VTK output file generation: material point data
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   04/04/2020
% Description:
% Function to generate a VTK file containing the material point data.
%
%--------------------------------------------------------------------------
% MAKEVTKMP(mpC,sig,uvw,mpFileName)
%--------------------------------------------------------------------------
% Input(s):
% mpC        - material point coordinates (nmp,nD)
% sig        - material point stresses (nmp,6)
% uvw        - material point displacements (nmp,nD)
% mpFileName - VTK file name, for example 'mpData.vtk'  
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

[nmp,nD]=size(mpC);                                                         % number of material points and dimensions

fid=fopen(mpFileName,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'MATLAB generated vtk file, WMC\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nmp);

%% position output 
if nD<3
    mpC = [mpC zeros(nmp,3-nD)];  
end
fprintf(fid,'%f %f %f \n',mpC');
fprintf(fid,'\n');

fprintf(fid,'POINT_DATA %i\n',nmp);

%% stress output
fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,1));
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,2));
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zz FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,3));
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_xy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,4));
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yz FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,5));
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zx FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,6));
fprintf(fid,'\n');

fprintf(fid,'SCALARS l_x FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',lp(:,1));
fprintf(fid,'\n');

fprintf(fid,'SCALARS l_y FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',lp(:,2));
fprintf(fid,'\n');

fprintf(fid,'SCALARS plastic FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',plastic(:));
fprintf(fid,'\n');


%% displacement output
if nD==3
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,2));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_z FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,3));
    fprintf(fid,'\n');
elseif nD==2
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,2));
    fprintf(fid,'\n');
elseif nD==1
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw);
    fprintf(fid,'\n');
end
fclose('all');
end