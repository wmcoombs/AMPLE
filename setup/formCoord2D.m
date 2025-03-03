function [etpl,coord,ftpl,fntpl,efnum] = formCoord2D(nelsx,nelsy,lx,ly)

%Two dimensional finite element grid generation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   06/05/2015
% Description:
% Function to generate a 2D finite element grid of linear quadrilateral 
% elements.
%
%--------------------------------------------------------------------------
% [etpl,coord,ftpl,fntpl,efnum] = FORMCOORD2D(nelsx,nelsy,lx,ly)
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
% ftpl  - face-element connections
% fntpl - face nodal topology
% efnum - local element face numbers associated with each face
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

if nargout>2
%% face generation (interior only)

nfaces = (nelsx-1)*nelsy + (nelsy-1)*nelsx;
ftpl   = zeros(nfaces,2);
fntpl  = zeros(nfaces,2);
efnum  = zeros(nfaces,2);
face   = 0;
for nely=1:nelsy
  for nelx=1:nelsx-1
    face=face+1;
    e1 = (nely-1)*(nelsx)+nelx;
    e2 = e1+1;
    ftpl(face,:) = [e1 e2];
    
    n1 = (nely-1)*(nelsx+1)+nelx+1;
    n2 = (nely)*(nelsx+1)+nelx+1;
    fntpl(face,:) = [n1 n2];
    
    efnum(face,:) = [1 3];
  end
end

for nelx=1:nelsx
  for nely=1:nelsy-1
    face=face+1;
    e1 = (nely-1)*(nelsx)+nelx;
    e2 = e1+nelsx;
    ftpl(face,:) = [e1 e2];
    
    n1 = (nely)*(nelsx+1)+nelx;
    n2 = n1+1;
    fntpl(face,:) = [n1 n2];    
    
    efnum(face,:) = [2 4];
  end
end

% figure;
% for i = 1:nfaces
%     plot(coord(fntpl(i,:)',1),coord(fntpl(i,:)',2),'b-'); hold on;
% end
% axis equal; axis off;

end
end