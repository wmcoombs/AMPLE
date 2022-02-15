function [nodes] = nodesForMP(etpl,elems)

%Unique list of nodes associated with a material point
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/05/2015
% Description:
% Function to determine a unique list of nodes for a group of elements
% connected to a material point.
%
%--------------------------------------------------------------------------
% [nodes] = NODESFORMP(etpl,elems)
%--------------------------------------------------------------------------
% Input(s):
% etpl  - element topology (nels,nen)
% elems - elements in the group (n,1)
%--------------------------------------------------------------------------
% Ouput(s);
% nodes - vector containing the nodes associated with the elements
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

nen=size(etpl,2);                                                           % number of nodes per element
n=size(elems,1)*size(elems,2);                                              % number of elements in group
nn=n*nen;                                                                   % number of nodes (inc. duplicates)
e=sort(reshape(etpl(elems,:),nn,1),1);                                      % list of all nodes (inc. duplicates)
d=[1; e(2:nn)-e(1:nn-1)]>0;                                                 % duplicate removal
nodes=e(d);                                                                 % unique list of nodes
end