function [lstps,g,mpData,mesh] = setupGrid_collapse

%Problem setup information
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Problem setup for a large deformation elasto-plastic collapse.
%
%--------------------------------------------------------------------------
% [lstps,g,mpData,mesh] = SETUPGRID2D
%--------------------------------------------------------------------------
% Input(s):
% 
%--------------------------------------------------------------------------
% Ouput(s);
% lstps  - number of loadsteps (1)
% g      - gravity (1)
% mpData - structured array with the following fields:
%           - mpType : material point type (1 = MPM, 2 = GIMPM)
%           - cmType : constitutive model type (1 = elastic, 2 = vM)
%           - mpC    : material point coordinates
%           - vp     : material point volume
%           - vp0    : initial material point volume
%           - mpM    : material point mass
%           - nIN    : nodes linked to the material point
%           - eIN    : element associated with the material point
%           - Svp    : basis functions for the material point
%           - dSvp   : basis function derivatives (at start of lstp)
%           - Fn     : previous deformation gradient
%           - F      : deformation gradient
%           - sig    : Cauchy stress
%           - epsEn  : previous logarithmic elastic strain
%           - epsE   : logarithmic elastic strain
%           - mCst   : material constants (or internal state parameters)
%           - fp     : force at the material point
%           - u      : material point displacement
%           - lp     : material point domain lengths
%           - lp0    : initial material point domain lengths
%
% mesh   - structured array with the following fields
%           - coord : mesh nodal coordinates (nodes,nD)
%           - etpl  : element topology (nels,nen)
%           - bc    : boundary conditions (*,2)
%           - h     : background mesh size (nD,1)
%--------------------------------------------------------------------------
% See also:
% FORMCOORD2D - background mesh generation
% DETMPPOS    - local material point positions
% SHAPEFUNC   - background grid basis functions
%--------------------------------------------------------------------------

%% Analysis parameters
E=1e6;   v=0.3;   fc=20e3;                                                  % Young's modulus, Poisson's ratio, yield strength   
mCst=[E v fc];                                                              % material constants
g=10;                                                                       % gravity
rho=1000;                                                                   % material density
lstps=40;                                                                   % number of loadsteps
a = 2;                                                                      % element multiplier
nelsx=4*a;                                                                  % number of elements in the x direction
nelsy=4*a;                                                                  % number of elements in the y direction
ly=8;  lx=8;                                                                % domain dimensions
mp=6;                                                                       % number of material points in each direction per element
mpType = 2;                                                                 % material point type: 1 = MPM, 2 = GIMP
cmType = 2;                                                                 % constitutive model: 1 = elastic, 2 = vM plasticity

%% Mesh generation
[etpl,coord] = formCoord2D(2*nelsx,nelsy,2*lx,ly);                          % background mesh generation
[nels,nen]   = size(etpl);                                                  % number of elements and nodes per element
[nodes,nD]   = size(coord);                                                 % number of nodes and dimensions
h            = [lx ly]./[nelsx nelsy];                                      % element lengths in each direction

%% Boundary conditions on backgroun mesh
bc = zeros(nodes*nD,2);                                                     % generate empty bc matrix
for node=1:nodes                                                            % loop over nodes
  if coord(node,1)==0                                                       % roller (x=0)
    bc(node*2-1,:)=[node*2-1 0];    
  end
  if coord(node,2)==0                                                       % roller (y=0)
    bc(node*2  ,:)=[node*2   0];
  end
end
bc = bc(bc(:,1)>0,:);                                                       % remove empty part of bc

%% Element limits for MP-element searching
eMin = zeros(nels,nD);                                                      % element lower coordinate limit
eMax = zeros(nels,nD);                                                      % element upper coordainte limit 
for i = 1:nD
    ci = coord(:,i);                                                        % nodal coordinates in current i direction
    c  = ci(etpl);                                                          % reshaped element coordinates in current i direction
    eMin(:,i) = min(c,[],2);                                                % element lower coordinate limit 
    eMax(:,i) = max(c,[],2);                                                % element upper coordainte limit 
end

%% Mesh data structure generation
mesh.etpl  = etpl;                                                          % element topology
mesh.coord = coord;                                                         % nodal coordinates
mesh.bc    = bc;                                                            % boundary conditions
mesh.h     = h;                                                             % mesh size
mesh.eMin  = eMin;                                                          % element lower coordinate limit 
mesh.eMax  = eMax;                                                          % element upper coordainte limit 

%% Material point generation
ngp    = mp^nD;                                                             % number of material points per element
GpLoc  = detMpPos(mp,nD);                                                   % local MP locations
N      = shapefunc(nen,GpLoc,nD);                                           % basis functions for the material points
[etplmp,coordmp] = formCoord2D(nelsx,nelsy,lx,ly);                          % mesh for MP generation
nelsmp = size(etplmp,1);                                                    % no. elements populated with material points
nmp    = ngp*nelsmp;                                                        % total number of mterial points

mpC=zeros(nmp,nD);                                                          % zero MP coordinates
for nel=1:nelsmp
  indx=(nel-1)*ngp+1:nel*ngp;                                               % MP locations within mpC
  eC=coordmp(etplmp(nel,:),:);                                              % element coordinates
  mpPos=N*eC;                                                               % global MP coordinates
  mpC(indx,:)=mpPos;                                                        % store MP positions
end
lp = zeros(nmp,2);                                                          % zero domain lengths
lp(:,1) = h(1)/(2*mp);                                                      % domain half length x-direction
lp(:,2) = h(2)/(2*mp);                                                      % domain half length y-direction
vp      = 2^nD*lp(:,1).*lp(:,2);                                            % volume associated with each material point

%% Material point structure generation
for mp = nmp:-1:1                                                           % loop backwards over MPs so array doesn't change size
  mpData(mp).mpType = mpType;                                               % material point type: 1 = MPM, 2 = GIMP
  mpData(mp).cmType = cmType;                                               % constitutive model: 1 = elastic, 2 = vM plasticity
  mpData(mp).mpC    = mpC(mp,:);                                            % material point coordinates
  mpData(mp).vp     = vp(mp);                                               % material point volume
  mpData(mp).vp0    = vp(mp);                                               % material point initial volume
  mpData(mp).mpM    = vp(mp)*rho;                                           % material point mass
  mpData(mp).nIN    = zeros(nen,1);                                         % nodes associated with the material point
  mpData(mp).eIN    = 0;                                                    % element associated with the material point
  mpData(mp).Svp    = zeros(1,nen);                                         % material point basis functions
  mpData(mp).dSvp   = zeros(nD,nen);                                        % derivative of the basis functions
  mpData(mp).Fn     = eye(3);                                               % previous deformation gradient
  mpData(mp).F      = eye(3);                                               % deformation gradient
  mpData(mp).sig    = zeros(6,1);                                           % Cauchy stress
  mpData(mp).epsEn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsE   = zeros(6,1);                                           % elastic strain (logarithmic)
  mpData(mp).mCst   = mCst;                                                 % material constants (or internal variables) for constitutive model
  mpData(mp).fp     = zeros(nD,1);                                          % point forces at material points
  mpData(mp).u      = zeros(nD,1);                                          % material point displacements
  if mpData(mp).mpType == 2
    mpData(mp).lp     = lp(mp,:);                                           % material point domain lengths (GIMP)
    mpData(mp).lp0    = lp(mp,:);                                           % initial material point domain lengths (GIMP)
  else
    mpData(mp).lp     = zeros(1,nD);                                        % material point domain lengths (MPM)
    mpData(mp).lp0    = zeros(1,nD);                                        % initial material point domain lengths (MPM)
  end
end
end
