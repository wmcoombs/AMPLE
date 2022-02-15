function makeVtk(coord,etpl,uvw,meshName)

%VTK output file generation: mesh data
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   03/04/2020
% Description:
% Function to generate a VTK file containing the background mesh data.
%
%--------------------------------------------------------------------------
% MAKEVTK(coord,etpl,meshName)
%--------------------------------------------------------------------------
% Input(s):
% coord    - coordinates of the grid nodes (nodes,nD)
% etpl     - element topology (nels,nen) 
% uvw      - incremental deformations (nodes*nD,1)
% meshName - VTK file name, for example 'mesh.vtk'  
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

[nodes,nD]=size(coord);
[nels,nen]=size(etpl);

%% FEM etpl to VTK format
if nD ==3
    if nen==20
        tvtk=[1 7 19 13 3 5 17 15 8 12 20 9 4 11 16 10 2 6 18 14];
        elemId=25;
        elemFormat='%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n';
    elseif nen==8
        tvtk=[1 4 8 5 2 3 7 6];
        elemId=12;
        elemFormat='%i %i %i %i %i %i %i %i %i\n';
    elseif nen==10
        tvtk=[1 2 3 4 5 6 7 8 10 9];
        elemId=24;
        elemFormat='%i %i %i %i %i %i %i %i %i %i %i\n';
    elseif nen==4
        tvtk=[1 3 2 4];
        elemId=10;
        elemFormat='%i %i %i %i %i\n';
    elseif nen==9
        tvtk=[3 1 7 5 2 8 6 4 9];
        elemId=10;
        elemFormat='%i %i %i %i %i %i %i %i %i %i\n';
    end
elseif nD==2
    if nen==3
        tvtk=[1 3 2];
        elemId=5;
        elemFormat='%i %i %i %i\n';
    elseif nen==4
        tvtk=[1 4 2 3];
        elemId=8;
        elemFormat='%i %i %i %i %i\n';
    elseif nen==8
        tvtk=[1 7 5 3 8 6 4 2];
        elemId=23;
        elemFormat='%i %i %i %i %i %i %i %i %i\n';
    end
end

%% Generation of vtk file
fid=fopen(meshName,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'MATLAB generated vtk file, WMC\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nodes);

%% nodal coordinates
if nD<3
    coord = [coord zeros(nodes,3-nD)];  
end
fprintf(fid,'%f %f %f \n',coord');
fprintf(fid,'\n');

%% element topology
fprintf(fid,'CELLS %i %i\n',nels,(nen+1)*nels);
etplOutput = [nen*ones(nels,1) (etpl(:,tvtk)-1)];
fprintf(fid,elemFormat,etplOutput');       
fprintf(fid,'\n');

%% element types
fprintf(fid,'CELL_TYPES %i\n',nels);
fprintf(fid,'%i\n',elemId*ones(nels,1));       
fprintf(fid,'\n');

%% displacement output
fprintf(fid,'POINT_DATA %i\n',nodes);
if nD==3
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(1:nD:end));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(2:nD:end));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_z FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(3:nD:end));
    fprintf(fid,'\n');
elseif nD==2
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(1:nD:end));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(2:nD:end));
    fprintf(fid,'\n');
elseif nD==1
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw);
    fprintf(fid,'\n');
end
fclose('all');
end