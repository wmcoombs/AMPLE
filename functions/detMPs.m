function [fint,Kt,mpData] = detMPs(uvw,mpData)

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
% FORMULSTIFF  - updated Lagrangian material stiffness calculation
% HOOKE3D      - linear elastic constitutive model
% VMCONST      - von Mises elasto-plastic constitutive model
%--------------------------------------------------------------------------

nmp   = length(mpData);
fint  = zeros(size(uvw));
npCnt = 0;                                                                  % counter for the number of entries in Kt
tnSMe = sum([mpData.nSMe]);                                                 % total number of stiffness matrix entries
krow  = zeros(tnSMe,1); kcol=krow; kval=krow;                               % zero the stiffness information

for mp=1:nmp                                                                % material point loop
    
    nIN     = mpData(mp).nIN;                                               % nodes associated with the material point 
    dNx     = mpData(mp).dSvp;                                              % basis function derivatives (start of lstp)
    [nD,nn] = size(dNx);                                                    % no. dimensions & no. nodes
    
    ed = reshape(ones(nD,1)*(nIN-1)*nD+(1:nD).'*ones(1,nn),1,nn*nD);        % degrees of freedom of nodes
    
    if nD==1                                                                % 1D case
        G=dNx;
        fPos=1;                                                             % deformation gradient positions
        aPos=1;                                                             % material stiffness matrix positions for global stiffness
        sPos=1;                                                             % Cauchy stress components for internal force
    elseif nD==2                                                            % 2D case (plane strain & stress)
        G=zeros(4,nD*nn);                                                   % zero the strain-disp matrix (2D)
        G([1 3],1:nD:end)=dNx;                                              % strain-displacement matrix
        G([4 2],2:nD:end)=dNx;
        fPos=[1 5 4 2];
        aPos=[1 2 4 5];
        sPos=[1 2 4 4];
    else                                                                    % 3D case
        G=zeros(9,nD*nn);                                                   % zero the strain-disp matrix (3D)
        G([1 4 9],1:nD:end)=dNx;
        G([5 2 6],2:nD:end)=dNx;
        G([8 7 3],3:nD:end)=dNx;
        fPos=[1 5 9 4 2 8 6 3 7];
        aPos=(1:9);
        sPos=[1 2 3 4 4 5 5 6 6];
    end
    
    ddF    = zeros(3); ddF(fPos) = G*uvw(ed);                               % spatial gradient (start of lstp) of displacements
    dF     = (eye(3)+ddF);                                                  % deformation gradient increment
    F      = dF*mpData(mp).Fn;                                              % deformation gradient
    epsEn  = mpData(mp).epsEn; epsEn(4:6) = 0.5*epsEn(4:6);                 % previous elastic strain
    epsEn  = epsEn([1 4 6; 4 2 5; 6 5 3]);                                  % matrix form of previous elastic strain
    [V,D]  = eig(epsEn);                                                    % eigen values and vectors of the elastic strain
    BeT    = dF*(V*diag(exp(2*diag(D)))*V.')*dF.';                          % trial left Cauchy-Green strain
    [V,D]  = eig(BeT);                                                      % eigen values and vectors of the trial left Cauchy-Green strain
    epsEtr = 0.5*V*diag(log(diag(D)))*V.';                                  % trial elastic strain (tensor form)
    epsEtr = [epsEtr([1 5 9]) 2*epsEtr([2 6 3])].';                         % trial elastic strain (vector form)
    
    %---------------------------------------------------------------------- % Constitutive model
    if mpData(mp).cmType == 1
        [D,Ksig,epsE]=Hooke3d(epsEtr,mpData(mp).mCst);                      % elastic behaviour
    elseif mpData(mp).cmType == 2
        [D,Ksig,epsE]=VMconst(epsEtr,mpData(mp).mCst);                      % elasto-plastic behaviour (von Mises)
    end
    %----------------------------------------------------------------------
    
    sig = Ksig/det(F);                                                      % Cauchy stress
    A   = formULstiff(F,D,sig,BeT);                                         % material stiffness
    
    dNx = dF(1:nD,1:nD)'\dNx;                                               % derivatives of basis functions (current)
    if nD==1                                                                % 1D case
        G=dNx;                                                              % strain-displacement matrix
    elseif nD==2                                                            % 2D case (plane strain & stress)
        G([1 3],1:nD:end)=dNx;                                              % strain-displacement matrix
        G([4 2],2:nD:end)=dNx;
    elseif nD==3                                                            % 3D case
        G([1 4 9],1:nD:end)=dNx;                                            % strain-displacement matrix
        G([5 2 6],2:nD:end)=dNx;
        G([8 7 3],3:nD:end)=dNx;
    end
    
    kp = mpData(mp).vp*det(dF)*(G.'*A(aPos,aPos)*G);                        % material point stiffness contribution
    fp = mpData(mp).vp*det(dF)*(G.'*sig(sPos));                             % internal force contribution
    
    mpData(mp).F    = F;                                                    % store deformation gradient
    mpData(mp).sig  = sig;                                                  % store Cauchy stress
    mpData(mp).epsE = epsE;                                                 % store elastic logarithmic strain
    
    npDoF=(size(ed,1)*size(ed,2))^2;                                        % no. entries in kp
    nnDoF=size(ed,1)*size(ed,2);                                            % no. DoF in kp
    krow(npCnt+1:npCnt+npDoF)=reshape(ed.'*ones(1,nnDoF),npDoF,1);          % row-column position storage
    kcol(npCnt+1:npCnt+npDoF)=reshape(ones(nnDoF,1)*ed  ,npDoF,1);
    kval(npCnt+1:npCnt+npDoF)=reshape(kp,npDoF,1);                          % stiffness storage
    npCnt=npCnt+npDoF;                                                      % number of entries in Kt
    fint(ed)=fint(ed)+fp;                                                   % internal force contribution
end

nDoF=length(uvw);                                                           % number of degrees of freedom
Kt=sparse(krow,kcol,kval,nDoF,nDoF);                                        % form the global stiffness matrix