function makeVtkMP(mpC,sig,uvw,mpFileName)

%VTK output file generation: material point data
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/01/2019
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

[nmp,nD]=size(mpC);

fid=fopen(mpFileName,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'MATLAB generated vtk file, WMC\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nmp);
if nD==3
    for i=1:nmp
        fprintf(fid,'%f %f %f \n',mpC(i,:));
    end
elseif nD==2
    for i=1:nmp
        fprintf(fid,'%f %f %f\n',mpC(i,:),0);
    end
elseif nD==1
    for i=1:nmp
        fprintf(fid,'%f %f %f\n',mpC(i),0,0);
    end
end
fprintf(fid,'\n');

fprintf(fid,'POINT_DATA %i\n',nmp);

fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',sig(i,1));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',sig(i,2));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zz FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',sig(i,3));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_xy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',sig(i,4));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yz FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',sig(i,5));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zx FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',sig(i,6));
end
fprintf(fid,'\n');

if nD==3
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',uvw(i,1));
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',uvw(i,2));
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_z FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',uvw(i,3));
    end
    fprintf(fid,'\n');
elseif nD==2
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',uvw(i,1));
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',uvw(i,2));
    end
    fprintf(fid,'\n');
elseif nD==1
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',uvw(i));
    end
    fprintf(fid,'\n');
end
fclose('all');