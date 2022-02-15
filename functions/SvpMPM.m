function [Svp,dSvp] = SvpMPM(xp,xv,h)

%1D material point basis functions
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   09/02/2016
% Description:
% Function to determine the one dimensional MPM shape functions based on
% global coordinates.
%
%--------------------------------------------------------------------------
% [Svp,dSvp] = SVPMPM(xp,xv,L)
%--------------------------------------------------------------------------
% Input(s):
% xp    - particle position
% xv    - grid node position
% h     - element length
%--------------------------------------------------------------------------
% Ouput(s);
% Svp   - particle characteristic function
% dSvp  - gradient of the characterstic function 
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

if     (xp-xv)<=-h                                                          % MP outside of element
  Svp =0; 
  dSvp=0;
elseif -h<(xp-xv) && (xp-xv)<=0                                             % MP in "left" element
  Svp = 1+(xp-xv)/h;
  dSvp= 1/h;
elseif  0<(xp-xv) && (xp-xv)<=h                                             % MP in "right" element
  Svp = 1-(xp-xv)/h;  
  dSvp=-1/h;
else                                                                        % MP outside of element
  Svp =0; 
  dSvp=0;
end