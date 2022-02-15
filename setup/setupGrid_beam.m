function [lstps,g,mpData,mesh] = setupGrid_beam

%Problem setup information
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Problem setup for a large deformation elastic beam analysis.
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
E=12e6;   v=0.2;   fc=20e4;                                                 % Young's modulus, Poisson's ratio, yield strength   
mCst=[E v fc];                                                              % material constants
g=10;                                                                       % gravity
rho=0;                                                                      % material density
P=-100e3;                                                                   % applied end load
lstps=50;                                                                   % number of loadsteps
a = 1;                                                                      % element multiplier
nelsx=22*a;                                                                 % number of elements in the x direction
nelsy=20*a;                                                                 % number of elements in the y direction
ly=10;  lx=11;                                                              % domain dimensions
d=1;  l=10;                                                                 % beam dimensions
mp=6;                                                                       % number of material points in each direction per element
mpType = 2;                                                                 % material point type: 1 = MPM, 2 = GIMP
cmType = 1;                                                                 % constitutive model: 1 = elastic, 2 = vM plasticity

%% Mesh generation
[etpl,coord] = formCoord2D(nelsx,nelsy,lx,ly);                              % background mesh generation
[~,nen]      = size(etpl);                                                  % number of element nodes
[nodes,nD]   = size(coord);                                                 % number of nodes and dimensions
h            = [lx ly]./[nelsx nelsy];                                      % element lengths in each direction

%% Boundary conditions on background mesh
bc = zeros(nodes*nD,2);                                                     % generate empty bc matrix
for node=1:nodes                                                            % loop over nodes
  if coord(node,1)==0                                                       % roller (x=0)
    bc(node*2-1,:)=[node*2-1 0];    
  end
  if coord(node,1)==0 && coord(node,2)==(ly-d/2)                            % mid-depth pin
    bc(node*2  ,:)=[node*2   0];
  end
end
bc = bc(bc(:,1)>0,:);                                                       % remove empty part of bc

%% Mesh data structure generation
mesh.etpl  = etpl;                                                          % element topology
mesh.coord = coord;                                                         % nodal coordinates
mesh.bc    = bc;                                                            % boundary conditions
mesh.h     = h;                                                             % mesh size

%% Material point generation
ngp    = mp^nD;                                                             % number of material points per element
GpLoc  = detMpPos(mp,nD);                                                   % local MP locations
N      = shapefunc(nen,GpLoc,nD);                                         % basis functions for the material points
[etplmp,coordmp] = formCoord2D(20*a,2*a,l,d);                               % mesh for MP generation
coordmp(:,2)=coordmp(:,2)+(ly-d);                                           % adjust MP locations (vertical)
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
  if abs(mpC(mp,1)-(l-lp(mp,1)))<1e-3/a && abs(mpC(mp,2)-(ly-d/2))<(lp(mp,2)+1e-3/a)
    mpData(mp).fp   = [0 P/2].';                                            % point forces at material points (end load)  
  else
    mpData(mp).fp   = zeros(nD,1);                                          % point forces at material points
  end
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