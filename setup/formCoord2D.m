function [etpl,coord] = formCoord2D(nelsx,nelsy,lx,ly)

%Two dimensional finite element grid generation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   06/05/2015
% Description:
% Function to generate a 2D finite element grid of linear quadrilateral 
% elements.
%
%--------------------------------------------------------------------------
% [etpl,coord] = FORMCOORD2D(nelsx,nelsy,lx,ly)
%--------------------------------------------------------------------------
% Input(s):
% nelsx - number of elements in the x direction
% nelsy - number of elements in the y direction
% lx    - length in the x direction
% ly    - length in the y direction
%--------------------------------------------------------------------------
% Ouput(s);
% etpl  - element topology
% coord - nodal coordinates
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

nels  = nelsx*nelsy;                                                        % number of elements
nodes = (nelsx+1)*(nelsy+1);                                                % number of nodes

%% node generation
coord = zeros(nodes,2);                                                     % zero coordinates 
node  = 0;                                                                  % zero node counter
for j=0:nelsy
  y=ly*j/nelsy;  
  for i=0:nelsx
    node=node+1;
    x=lx*i/nelsx;
    coord(node,:)=[x y];
  end
end
%% element generation
etpl = zeros(nels,4);                                                       % zero element topology
nel  = 0;                                                                   % zero element counter
for nely=1:nelsy
  for nelx=1:nelsx
    nel=nel+1;
    etpl(nel,1)=(nely-1)*(nelsx+1)+nelx;
    etpl(nel,2)=nely*(nelsx+1)+nelx;  
    etpl(nel,3)=etpl(nel,2)+1;
    etpl(nel,4)=etpl(nel,1)+1; 
  end
end
end