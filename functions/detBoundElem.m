function [mesh] = detBoundElem(mesh)

%Boundary element determination
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   11/06/2021
% Description:
% Function to determine the bounadry elements and the faces attached to the
% boundary elements. 
%
%--------------------------------------------------------------------------
% [mesh] = DETBOUNDELEM(mesh)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh structured array. Function requires: 
%           - coord : coordinates of the grid nodes (nodes,nD)
%           - etpl  : element topology (nels,nen) 
%           - ftpl  : face-element interactions
%--------------------------------------------------------------------------
% Ouput(s);
% mesh   - mesh structured array. Function modifies:
%           - bE    : bounadry elements flag
%           - bF    : boundary face flag 
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------


ftpl  = mesh.ftpl;                                                          % face topology
eInA  = mesh.eInA;                                                          % elements in analysis
nface = length(ftpl);                                                       % number of faces
A  = eInA(ftpl);                                                            % elements attached to faces in/out
B  = ne(A(:,1),A(:,2));                                                     % faces between in/out elements
bF = (1:nface).';                                                           % list, 1 to number of faces
bF = bF(B);                                                                 % face numbers of faces between in and out elements
bE = zeros(size(bF));                                                       % bounday element storage
for i = 1:length(bF)                                                        % loop over in/out faces 
    fE = ftpl(bF(i),:);                                                     % elements attached to face                                                 
    if eInA(fE(1))==1
        bE(i) = fE(1);                                                      % set boundary element = element 1
    else 
        bE(i) = fE(2);                                                      % set boundary element = element 2
    end
end
bE = bE(bE>0);
bE = unique(bE);                                                            % remove duplicates
a  = false(nface,1);                                                        % zero boundary face flag
for i = 1:nface
    fE = ftpl(i,:);                                                         % elements attached to face
    b     = sum(eInA(fE))==2;                                               % both elements in analysis (yes, b1=1; no, b=0)
    a(i)  = (sum(bE==fE(1))+sum(bE==fE(2)))*b>=1;                           % is face attached to boundary element 
end
mesh.bE = bE;                                                               % store boundary elements
mesh.bF = a;                                                                % store boundary face flag

end
