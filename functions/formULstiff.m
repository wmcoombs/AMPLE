function [A] = formULstiff(F,D,s,B)

%Updated Lagrangian material stiffness matrix
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   27/05/2015
% Description:
% Function to determine consistent material stiffness matrix based on an
% updated Lagrangian formulation of finite deformation mechanics.  See
% equations (25) and (26) of the following paper for full details:
% Charlton, T.J., Coombs, W.M. & Augarde, C.E. (2017). iGIMP: An implicit 
% generalised interpolation material point method for large deformations. 
% Computers and Structures 190: 108-125.
%
%--------------------------------------------------------------------------
% [A] = FORMULSTIFF(F,D,s,BeT)
%--------------------------------------------------------------------------
% Input(s):
% F  - deformation gradient (3,3)
% D  - small strain material stifness matrix (6,6)
% s  - Cauchy stress (6,1)
% B  - trial elastic left Cauchy-Green strain matrix (3,3)
%--------------------------------------------------------------------------
% Ouput(s);
% A   - consistent tangent stifness matrix (9,9)
%--------------------------------------------------------------------------
% See also:
% PARDERGEN  - partial derivative of a second order tensor
%--------------------------------------------------------------------------

t = [1:4 4 5 5 6 6];                                                        % 6 to 9 component steering vector
J = det(F);                                                                 % volume ratio
[bV,bP] = eig(B); bP = [bP(1); bP(5); bP(9)];                               % eigen values/vector of the trial elastic strain tensor
L = parDerGen(B,bV,bP,log(bP),1./bP);                                       % derivative of the logarithmic strain
S = zeros(9); T = S;                                                        % zero matrices
S(1,[1 4 9]) = s([1 4 6]);                                                  % matrix form of sigma_{il}delta_{jk}
S(2,[2 5 6]) = s([2 4 5]); 
S(3,[3 7 8]) = s([3 5 6]);  
S(4,[2 5 6]) = s([4 1 6]); 
S(5,[1 4 9]) = s([4 2 5]); 
S(6,[3 7 8]) = s([5 2 4]);
S(7,[2 5 6]) = s([5 6 3]); 
S(8,[1 4 9]) = s([6 5 3]); 
S(9,[3 7 8]) = s([6 4 1]);
T(1,[1 4 9])     = 2*B([1 4 7]);                                            % matrix form of delta_{pk}b^e_{ql}+delta_{qk}b^e_{pl}+
T(2,[2 5 6])     = 2*B([5 2 8]); 
T(3,[3 7 8])     = 2*B([9 6 3]); 
T(4,[1 2 4:6 9]) = B([2 4 5 1 7 8]); 
T(5,:)           = T(4,:);
T(6,[2 3 5:8  ]) = B([6 8 3 9 5 2]); 
T(7,:)           = T(6,:);
T(8,[1 3 4 7:9]) = B([3 7 6 4 1 9]); 
T(9,:)           = T(8,:);
A = (D(t,t))*L(t,t)*T/(2*J)-S;                                              % consistent tangent stiffness matrix




function [L] = parDerGen(X,eV,eP,yP,ydash)

%Partial derivative of a second order tensor function
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   27/05/2015
% Description:
% Function to determine the partial derivative of a second order tensor
% function with respect to its arguement (X) based on the implementation 
% described by in the following paper:
%
% C. Miehe, Comparison of two algorithms for the computation of fourth-
% order isotropic tensor functions, Computers & Structures 66 (1998) 37-43.
%
% For example, in order to determine the derivative of log(X) with respect
% to X the inputs to the function should be:
%
% [L] = PARDERGEN(X,eV,eP,log(eP),1./eP)
%
% as the derivative of the log(x) is 1/x
%
% The symbols used in the code follow, as closely as possible, those used
% in the Miehe (1998) paper.  There are a number of different cases that
% have to be checked (all zero and repeated eigenvalues) in addition to the
% general case where there are no repeated eigenvalues.  
%
%--------------------------------------------------------------------------
% [L] = PARDERGEN(X,eV,eP,yP,ydash)
%--------------------------------------------------------------------------
% Input(s):
% X     - second order tensor in matrix format (3,3)
% eV    - eigenvectors of X (3,3)
% eP    - eigenvalues of X (1,3) 
% yP    - function applied to eP (1,3)
% ydash - derivative of the function applied to eP (1,3)
%--------------------------------------------------------------------------
% Ouput(s);
% L      - partial derivative of the second order tensor with respect to 
%          its arguement (6,6)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------


tol=1e-9;
Is=[eye(3) zeros(3); zeros(3) eye(3)/2];
bm1=[1 1 1 0 0 0].';
if (abs(eP(1))<tol && abs(eP(2))<tol && abs(eP(3))<tol)                     % all zero eigenvalues case
    L = Is;
elseif abs(eP(1)-eP(2))<tol && abs(eP(1)-eP(3))<tol                         % equal eigenvalues case
    L = ydash(1)*Is;
elseif abs(eP(1)-eP(2))<tol || abs(eP(2)-eP(3))<tol || abs(eP(1)-eP(3))<tol % repeated eigenvalues case
    if     abs(eP(1)-eP(2))<tol
        xa  = eP(3);    xc  = eP(1);
        ya  = yP(3);    yc  = yP(1);
        yda = ydash(3); ydc = ydash(1);
    elseif abs(eP(2)-eP(3))<tol
        xa  = eP(1);    xc  = eP(2);
        ya  = yP(1);    yc  = yP(2);
        yda = ydash(1); ydc = ydash(2);
    else
        xa  = eP(2);    xc  = eP(1);
        ya  = yP(2);    yc  = yP(1);
        yda = ydash(2); ydc = ydash(1);
    end
    x    = X([1 5 9 4 6 3]).';
    s    = zeros(5,1);
    s(1) = (ya-yc)/(xa-xc)^2-ydc/(xa-xc);
    s(2) = 2*xc*(ya-yc)/(xa-xc)^2-(xa+xc)/(xa-xc)*ydc;
    s(3) = 2*(ya-yc)/(xa-xc)^3-(yda+ydc)/(xa-xc)^2;
    s(4) = xc*s(3);
    s(5) = xc^2*s(3);
    dX2dX=[2*X(1) 0      0      X(2)         0             X(3)          ;
           0      2*X(5) 0      X(2)         X(6)          0             ;
           0      0      2*X(9) 0            X(6)          X(3)          ;
           X(2)   X(2)   0     (X(1)+X(5))/2 X(3)/2        X(6)/2        ;
           0      X(6)   X(6)   X(3)/2       (X(5)+X(9))/2 X(2)/2        ;
           X(3)   0      X(3)   X(6)/2       X(2)/2        (X(1)+X(9))/2];
    L = s(1)*dX2dX-s(2)*Is-s(3)*(x*x.')+s(4)*(x*bm1.'+bm1*x.')-s(5)*(bm1*bm1.');
else                                                                        % general case (non repeated eigenvalues)
    D=[(eP(1)-eP(2))*(eP(1)-eP(3));
       (eP(2)-eP(1))*(eP(2)-eP(3));
       (eP(3)-eP(1))*(eP(3)-eP(2))];
    alfa=0; bta=0; gama=zeros(3,1); eDir=zeros(6,3);
    for i=1:3
        alfa = alfa+yP(i)*eP(i)/D(i);
        bta  = bta+yP(i)/D(i)*det(X);
        for j=1:3
            gama(i) = gama(i)+yP(j)*eP(j)/D(j)*(det(X)/eP(j)-eP(i)^2)*1/eP(i)^2;
        end
        esq = eV(:,i)*eV(:,i).';
        eDir(:,i) = [esq(1,1) esq(2,2) esq(3,3) esq(1,2) esq(2,3) esq(3,1)].';
    end
    y = inv(X);
    Ib=[y(1)^2    y(2)^2    y(7)^2     y(1)*y(2)               y(2)*y(7)               y(1)*y(7)              ;
        y(2)^2    y(5)^2    y(6)^2     y(5)*y(2)               y(5)*y(6)               y(2)*y(6)              ;
        y(7)^2    y(6)^2    y(9)^2     y(6)*y(7)               y(9)*y(6)               y(9)*y(7)              ;
        y(1)*y(2) y(5)*y(2) y(6)*y(7) (y(1)*y(5)+y(2)^2)/2    (y(2)*y(6)+y(5)*y(7))/2 (y(1)*y(6)+y(2)*y(7))/2 ;
        y(2)*y(7) y(5)*y(6) y(9)*y(6) (y(2)*y(6)+y(5)*y(7))/2 (y(9)*y(5)+y(6)^2)/2    (y(9)*y(2)+y(6)*y(7))/2 ;
        y(1)*y(7) y(2)*y(6) y(9)*y(7) (y(1)*y(6)+y(2)*y(7))/2 (y(9)*y(2)+y(6)*y(7))/2 (y(9)*y(1)+y(7)^2)/2   ];
    L = alfa*Is-bta*Ib;
    for i=1:3
        L = L+(ydash(i)+gama(i))*eDir(:,i)*eDir(:,i).';
    end
end