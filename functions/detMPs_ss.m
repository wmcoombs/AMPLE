function [fint,Kt,mpData] = detMPs_ss(uvw,mpData)

%Stiffness and internal force calculation for all material points
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to determine the stiffness contribution of a particle to the
% nodes that it influences based on a Updated Lagrangian finite deformation 
% formulation.  The function also returns the stresses at the particles and 
% the internal force contribution.  This function allows for elasto-
% plasticity at the material points.  The functionis applicable to 1, 2 and
% 3 dimensional problems without modification as well as different material 
% point methods and background meshes.   
% 
%--------------------------------------------------------------------------
% [fint,Kt,mpData] = DETMPS(uvw,mpData)
%--------------------------------------------------------------------------
% Input(s):
% uvw    - nodal displacements that influence the MP (nn*nD,1)
% mpData - material point structured array. The following fields are
%          required by the function:
%           - dSvp  : basis function derivatives (nD,nn)
%           - nIN   : background mesh nodes associated with the MP (1,nn)
%           - Fn    : previous deformation gradient (3,3) 
%           - epsEn : previous elastic logarithmic strain (6,1)
%           - mCst  : material constants
%           - vp    : material point volume (1)
%           - nSMe  : number stiffness matrix entries
% nD     - number of dimensions
%--------------------------------------------------------------------------
% Ouput(s);
% fint   - global internal force vector
% Kt     - global stiffness matrix
% mpData - material point structured array (see above).  The following
%          fields are updated by the function:
%           - F     : current deformation gradient (3,3)
%           - sig   : current Cauchy stress (6,1)
%           - epsE  : current elastic logarithmic strain (6,1)
%--------------------------------------------------------------------------
% See also:
% FORMULSTIFF      - updated Lagrangian material stiffness calculation
% HOOKE3D          - linear elastic constitutive model
% VMCONST          - von Mises elasto-plastic constitutive model
%--------------------------------------------------------------------------

nmp   = length(mpData);                                                     % number of material points
fint  = zeros(size(uvw));                                                   % zero internal force vector
npCnt = 0;                                                                  % counter for the number of entries in Kt
tnSMe = sum([mpData.nSMe]);                                                 % total number of stiffness matrix entries
krow  = zeros(tnSMe,1); kcol=krow; kval=krow;                               % zero the stiffness information
ddF   = zeros(3);                                                           % derivative of duvw wrt. spatial position

nD = length(mpData(1).mpC);

if nD==1                                                                    % 1D case
    aPos=1;                                                                 % material stiffness matrix positions for global stiffness
    sPos=1;                                                                 % Cauchy stress components for internal force
elseif nD==2                                                                % 2D case (plane strain & stress)
    aPos=[1 2 4];
    sPos=[1 2 4];
else                                                                        % 3D case
    aPos=[1 2 3 4 5 6];
    sPos=[1 2 3 4 5 6];
end

for mp=1:nmp                                                                % material point loop
    
    nIN = mpData(mp).nIN;                                                   % nodes associated with the material point 
    dNx = mpData(mp).dSvp;                                                  % basis function derivatives (start of lstp)
    nn  = size(dNx,2);                                                      % no. dimensions & no. nodes
    ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                    % degrees of freedom of nodes (matrix form)
    ed  = reshape(ed,1,nn*nD);                                              % degrees of freedom of nodes (vector form)
    
    if nD==1                                                                % 1D case
        G=dNx;                                                              % strain-displacement matrix
    elseif nD==2                                                            % 2D case (plane strain & stress)
        G=zeros(3,nD*nn);                                                   % zero the strain-disp matrix (2D)
        G([1 3],1:nD:end)=dNx;                                              % strain-displacement matrix
        G([3 2],2:nD:end)=dNx;
    else                                                                    % 3D case
        G=zeros(6,nD*nn);                                                   % zero the strain-disp matrix (3D)
        G([1 4 6],1:nD:end)=dNx;                                            % strain-displacement matrix
        G([4 2 5],2:nD:end)=dNx;
        G([6 5 3],3:nD:end)=dNx;
    end
    
    deps = zeros(6,1);
    deps(sPos)   = G*uvw(ed);                                              %  strain increment
    epsEn  = mpData(mp).epsEn;                                              % previous elastic strain
    epsEtr = epsEn + deps;                                                  % trial elastic strain
    
    %---------------------------------------------------------------------- % Constitutive model
    if mpData(mp).cmType == 1
        [A,sig,epsE]=Hooke3d(epsEtr,mpData(mp).mCst);                      % elastic behaviour
    elseif mpData(mp).cmType == 2
        [A,sig,epsE]=VMconst(epsEtr,mpData(mp).mCst);                      % elasto-plastic behaviour (von Mises)
    end
    %----------------------------------------------------------------------
    
    kp = mpData(mp).vp*(G.'*A(aPos,aPos)*G);                                % material point stiffness contribution
    fp = mpData(mp).vp*(G.'*sig(sPos));                                     % internal force contribution
    
    mpData(mp).sig  = sig;                                                  % store Cauchy stress
    mpData(mp).epsE = epsE;                                                 % store elastic logarithmic strain
    
    npDoF=(size(ed,1)*size(ed,2))^2;                                        % no. entries in kp
    nnDoF=size(ed,1)*size(ed,2);                                            % no. DoF in kp                        
    krow(npCnt+1:npCnt+npDoF)=repmat(ed.',nnDoF,1);                         % row position storage
    kcol(npCnt+1:npCnt+npDoF)=repmat(ed  ,nnDoF,1);                         % column position storage
    kval(npCnt+1:npCnt+npDoF)=kp;                                           % stiffness storage
    npCnt=npCnt+npDoF;                                                      % number of entries in Kt
    fint(ed)=fint(ed)+fp;                                                   % internal force contribution
end

nDoF=length(uvw);                                                           % number of degrees of freedom
Kt=sparse(krow,kcol,kval,nDoF,nDoF);                                        % form the global stiffness matrix
end