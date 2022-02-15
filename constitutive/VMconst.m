function [Dalg,sig,epsE] = VMconst(epsEtr,mCst)

%von Mises linear elastic perfectly plastic constitutive model
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   16/05/2016
% Description:
% von Mises perfect plasticity constitutive model with an implicit backward
% Euler stress integration algorithm based on the following thesis:
%
% Coombs, W.M. (2011). Finite deformation of particulate geomaterials: 
% frictional and anisotropic Critical State elasto-plasticity. School of 
% Engineering and Computing Sciences. Durham University. PhD.
%
%--------------------------------------------------------------------------
% [Dalg,sigma,epsE] = VMCONST(epsEtr,mCst)
%--------------------------------------------------------------------------
% Input(s):
% epsEtr - trial elastic strain (6,1)
% mCst   - material constants 
%--------------------------------------------------------------------------
% Ouput(s);
% sig    - Cauchy stress (6,1)
% epsE   - elastic strain (6,1)
% Dalg   - algorithmic consistent tangent (6,6)
%--------------------------------------------------------------------------
% See also:
% YILEDFUNCDERIVATIVES - yield function 1st and 2nd derivatives
%--------------------------------------------------------------------------

E=mCst(1); v=mCst(2); rhoY=mCst(3);                                         % material constants
tol=1e-9; maxit=5 ;                                                         % NR parameters
bm1=[1 1 1 0 0 0]';                                                         % vector form of an identity matrix
Ce=[-ones(3)*v+(1+v)*eye(3) zeros(3); zeros(3) 2*(1+v)*eye(3)]/E;           % elastic compliance matrix
De=E/((1+v)*(1-2*v))*(bm1*bm1'*v+...
    [eye(3) zeros(3); zeros(3) eye(3)/2]*(1-2*v));                          % elastic stiffness matrix
sig=De*epsEtr;                                                              % elastic trial stress
s=sig-sum(sig(1:3))/3*bm1; j2=(s'*s+s(4:6)'*s(4:6))/2;                      % deviatoric stress and its second invariant (J2)
f=sqrt(2*j2)/rhoY-1;                                                        % yield function
epsE=epsEtr; Dalg=De;                                                       % set the elastic case
if (f>tol)                                                                  % plasticity loop
  b=zeros(7,1); b(7)=f; itnum=0; dgam=0;                                    % initial conditions on the NR search
  [df,ddf] = yieldFuncDerivatives(sig,rhoY);                                % 1st and 2nd derivative of f wrt. stress
  while (itnum<maxit) && ((norm(b(1:6))>tol) || (abs(b(7))>tol))            % NR loop
    A=[eye(6)+dgam*ddf*De df; df'*De 0];                                    % derivative of the residuals wrt. unknowns (Hessian)
    dx=-inv(A)*b; epsE=epsE+dx(1:6); dgam=dgam+dx(7);                       % increment unknowns
    sig=De*epsE;                                                            % updated stress
    s=sig-sum(sig(1:3))/3*bm1; j2=(s'*s+s(4:6)'*s(4:6))/2;                  % deviatoric stress and its second invariant (J2)  
    [df,ddf] = yieldFuncDerivatives(sig,rhoY);                              % 1st and 2nd derivative of f wrt. stress
    b=[epsE-epsEtr+dgam*df; sqrt(2*j2)/rhoY-1];                             % residuals
    itnum=itnum+1;                                                          % iteration counter
  end
  B=inv([Ce+dgam*ddf df; df.' 0]);                                          % Aalg from eqn (2.53) of Coombs(2011) thesis
  Dalg=B(1:6,1:6);                                                          % algorithmic consistent tangent
end
end

function [df,ddf] = yieldFuncDerivatives(sig,rhoY)

%von Mises yield function derivatives
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   16/05/2016
% Description:
% First and second derivatives of the von Mises yield function with respect
% to stress.
%
%--------------------------------------------------------------------------
% [df,ddf] = YIELDFUNCDERIVATIVES(sigma,rhoY)
%--------------------------------------------------------------------------
% Input(s):
% sigma - Cauchy stress (6,1)
% rhoY  - von Mises yield strength (1)
%--------------------------------------------------------------------------
% Ouput(s);
% df    - derivative of the yield function wrt. sigma (6,1)
% ddf   - second derivative of the yield function wrt. sigma (6,6)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

bm1=[1 1 1 0 0 0]';                                                         % vector form of an identity matrix
s=sig-sum(sig(1:3))/3*bm1; j2=(s'*s+s(4:6)'*s(4:6))/2;                      % deviatoric stress and its second invariant (J2)
dj2=s; dj2(4:6)=2*dj2(4:6);                                                 % derivative of J2 wrt. stress
ddj2=[eye(3)-ones(3)/3 zeros(3); zeros(3) 2*eye(3)];                        % 2nd derivative of J2 wrt. stress
df =dj2/(rhoY*sqrt(2*j2));                                                  % derivative of f wrt. stress
ddf=1/rhoY*(ddj2/sqrt(2*j2)-(dj2*dj2.')/(2*j2)^(3/2));                      % 2nd derivative of f wrt. stress
end