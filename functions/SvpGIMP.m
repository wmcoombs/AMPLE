function [Svp,dSvp] = SvpGIMP(xp,xv,h,lp)

%1D generalised interpolation material point basis functions
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/05/2015
% Description:
% Function to determine the one dimensional GIMP shape functions based on
% the paper: 
% S.G. Bardenhagen, E.M. Kober, The generalized interpolation material 
% point method, Computer Modeling in Eng. & Sciences 5 (2004) 477-496.
%
% The letters on the different cases refer to the regions show in the
% AMPLE paper.  
%--------------------------------------------------------------------------
% [Svp,dSvp] = SVPGIMP(xp,xv,h,lp)
%--------------------------------------------------------------------------
% Input(s):
% xp    - particle position
% xv    - grid node position
% h     - grid spacing 
% lp    - particle half width 
%--------------------------------------------------------------------------
% Ouput(s);
% Svp   - particle characteristic function
% dSvp  - gradient of the characterstic function 
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

if     (-h-lp)<(xp-xv) && (xp-xv)<=(-h+lp)                                  % A: partial overlap of domain and element
  Svp = (h+lp+(xp-xv))^2/(4*h*lp);  
  dSvp= (h+lp+(xp-xv))/(2*h*lp);
elseif (-h+lp)<(xp-xv) && (xp-xv)<=(  -lp)                                  % B: full overlap of domain and element (left)
  Svp = 1+(xp-xv)/h;
  dSvp= 1/h;
elseif (  -lp)<(xp-xv) && (xp-xv)<=(   lp)                                  % C: partial overlap of domain and two elements element
  Svp = 1-((xp-xv)^2+lp^2)/(2*h*lp);
  dSvp=-(xp-xv)/(h*lp);
elseif (   lp)<(xp-xv) && (xp-xv)<=( h-lp)                                  % D: full overlap of domain and element (right)
  Svp = 1-(xp-xv)/h;  
  dSvp=-1/h;
elseif ( h-lp)<(xp-xv) && (xp-xv)<=( h+lp)                                  % E: partial overlap of domain and element                    
  Svp = (h+lp-(xp-xv))^2/(4*h*lp); 
  dSvp=-(h+lp-(xp-xv))/(2*h*lp);
else                                                                        % zero overlap
  Svp =0; 
  dSvp=0;
end