function topo = get_adj_mat(V, C)
% GET_ADJ_MAT  Build the full topology (adjacency) description of the
% vertex model tissue from its vertices and cell-vertex lists.
%
%   topo = GET_ADJ_MAT(V, C) computes the vertex-vertex, vertex-edge,
%   cell-vertex, and cell-edge adjacency matrices, plus the edge list,
%   from the raw geometric description of the tissue, and packages them
%   into a single "topo" struct used throughout the rest of the code.
%
% Inputs:
%   V - Nvert x 2 array of vertex positions
%   C - Ncell x 1 cell array; C{i} lists the vertex indices (in order
%       around the polygon) belonging to cell i
%
% Output:
%   topo - struct with fields:
%       Ncell, Nvert, Nedge - counts of cells, vertices, edges
%       VVadj  - Nvert x Nvert logical vertex-vertex adjacency matrix
%       E_list - Nedge x 2 list of vertex-index pairs forming each edge
%       VEadj  - Nvert x Nedge vertex-edge incidence matrix
%       CVadj  - Ncell x Nvert cell-vertex incidence matrix
%       CEadj  - Ncell x Nedge cell-edge incidence matrix (logical)

Ncell = size(C, 1);
Nvert = size(V, 1);

% --- Vertex-vertex adjacency: connect consecutive vertices around each
% cell polygon (closing the loop from the last vertex back to the first).
VVadj = zeros(Nvert, Nvert);
for i = 1:Ncell
    idx = C{i};
    ZC = length(idx);
    for j = 1:(ZC-1)
        VVadj(idx(j), idx(j+1)) = 1;
    end
    VVadj(idx(ZC), idx(1)) = 1;
end
VVadj = VVadj + transpose(VVadj);   % symmetrize (undirected adjacency)
VVadj = logical(VVadj);

% --- Edge list: one row per unique vertex pair (upper triangle only, to
% avoid double-counting each undirected edge).
[vv1, vv2] = find(triu(VVadj));
E_list = [vv1, vv2];
Nedge = size(E_list, 1);

% --- Vertex-edge incidence: mark which two vertices bound each edge.
VEadj = zeros(Nvert, Nedge);
for eid = 1:Nedge
    VEadj(vv1(eid), eid) = 1;
    VEadj(vv2(eid), eid) = 1;
end

% --- Cell-vertex incidence: mark which vertices belong to each cell.
CVadj = zeros(Ncell, Nvert);
for i = 1:Ncell
    CVadj(i, C{i}) = 1;
end

% --- Cell-edge incidence: for each cell, find the edge connecting each
% pair of consecutive vertices around its polygon (closing the loop from
% the last vertex back to the first).
CEadj = zeros(Ncell, Nedge);
for i = 1:Ncell
    zz = length(C{i});
    for j = 1:(zz-1)
        eid = intersect(find(VEadj(C{i}(j), :)), find(VEadj(C{i}(j+1), :)));
        CEadj(i, eid) = 1;
    end
    eid = intersect(find(VEadj(C{i}(zz), :)), find(VEadj(C{i}(1), :)));
    CEadj(i, eid) = 1;
end
CEadj = logical(CEadj);

% --- Package everything into the topology struct used elsewhere.
topo = struct(...
    'Ncell', Ncell, ...
    'Nvert', Nvert, ...
    'Nedge', Nedge, ...
    'VVadj', VVadj, ...
    'CEadj', CEadj, ...
    'CVadj', CVadj, ...
    'E_list', E_list, ...
    'VEadj', VEadj ...
    );
end
