function [D,sig,epsE] = Hooke3d(epsE,mCst)
%Linear elastic constitutive model
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Small strain linear elastic constitutive model 
%
%--------------------------------------------------------------------------
% [D,sig,epsE] = HOOKE3D(epsE,mCst)
%--------------------------------------------------------------------------
% Input(s):
% epsE  - elastic strain (6,1)
% mCst  - material constants 
%--------------------------------------------------------------------------
% Ouput(s);
% D     - elastic stiffness matrix (6,6)
% sig   - stress (6,1)
% epsE  - elastic strain (6,1)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

E=mCst(1);                                                                  % Young's modulus
v=mCst(2);                                                                  % Poisson's ratio
bm1=[1 1 1 0 0 0]';
D=E/((1+v)*(1-2*v))*(bm1*bm1'*v+...                                         % elastic stiffness
  [eye(3) zeros(3); zeros(3) eye(3)/2]*(1-2*v)); 
sig=D*epsE;                                                                 % stress 