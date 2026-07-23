function [V, C, topo] = T1_edgeswap_LE(V, C, topo, T1_edge_ID, strain)
% T1_EDGESWAP_LE  Perform a T1 (neighbor-exchange) topological transition
% on a single edge, under Lees-Edwards (sheared) periodic boundary
% conditions.
%
%   [V, C, topo] = T1_EDGESWAP_LE(V, C, topo, T1_edge_ID, strain)
%   collapses the edge T1_edge_ID (connecting vertices u0 and v0, shared
%   by cells c1 and c2) and re-expands it perpendicular to its original
%   orientation, so that the two cells that used to touch along this edge
%   (c1, c2) separate, and the two cells that used to only meet at a
%   point (cu, cv) become new neighbors. Same transition as
%   T1_edgeswap.m, but the edge geometry is unwrapped/re-wrapped using
%   the sheared PBC convention.
%
% Inputs:
%   V          - Nvert x 2 array of vertex positions
%   C          - Ncell x 1 cell array of vertex indices per cell
%   topo       - topology struct (see get_adj_mat)
%   T1_edge_ID - index (row of topo.E_list) of the edge to flip
%   strain     - current Lees-Edwards shear strain
%
% Outputs:
%   V, C, topo - updated with the new vertex positions and topology
%                after the T1 transition

box_size = sqrt(topo.Ncell);

u0 = topo.E_list(T1_edge_ID, 1);
v0 = topo.E_list(T1_edge_ID, 2);

% --- Identify the two cells (c1, c2) sharing this edge, and orient them
% consistently: c1 is the cell for which going u0 -> v0 follows the
% cell's own (counter)clockwise vertex ordering.
temp = find(topo.CEadj(:, T1_edge_ID));
c1 = temp(1);
c2 = temp(2);
zz = length(C{c1});
ind1 = find(C{c1} == u0);
if C{c1}((1 + mod(ind1+1-1, zz))) == v0
    % already correctly oriented
else
    c1 = temp(2);
    c2 = temp(1);
end

z1 = length(C{c1});
z2 = length(C{c2});
if z1 < 4 || z2 < 4
    % A T1 would reduce one of these cells below a triangle (3 vertices),
    % which is not a valid polygon, so the transition is skipped.
    disp('there are triangles, no T1 possible!');
    return;
end

% --- Find the two "outside" cells cu, cv that meet at vertices u0, v0
% respectively but are not c1 or c2 -- these will gain the new edge.
cell_neibs = topo.CVadj(:, u0);
cell_neibs([c1, c2]) = 0;
cu = find(cell_neibs);

cell_neibs = topo.CVadj(:, v0);
cell_neibs([c1, c2]) = 0;
cv = find(cell_neibs);

% --- Update cell cu's vertex list: insert v0 next to u0 (cu gains a
% vertex, since it will now also touch v0).
cell_list_cu = C{cu};
zz = length(cell_list_cu);
ind0 = find(cell_list_cu == u0);
ind1 = (1 + mod(ind0-1-1, zz));
ind2 = (1 + mod(ind0+1-1, zz));
u1 = cell_list_cu(ind1); %#ok<NASGU> % neighbor vertices, kept for clarity
u2 = cell_list_cu(ind2);
cell_list_cu = [cell_list_cu(1:ind0), v0, cell_list_cu(ind0+1:end)];

% --- Update cell cv's vertex list: insert u0 next to v0 (symmetric to cu).
cell_list_cv = C{cv};
zz = length(cell_list_cv);
ind0 = find(cell_list_cv == v0);
ind1 = (1 + mod(ind0-1-1, zz));
ind2 = (1 + mod(ind0+1-1, zz));
v1 = cell_list_cv(ind1); %#ok<NASGU>
v2 = cell_list_cv(ind2);
cell_list_cv = [cell_list_cv(1:ind0), u0, cell_list_cv(ind0+1:end)];

% --- Cells c1, c2 each lose one vertex: c1 no longer touches u0, c2 no
% longer touches v0.
cell_list_c1 = C{c1};
cell_list_c1(cell_list_c1 == u0) = [];

cell_list_c2 = C{c2};
cell_list_c2(cell_list_c2 == v0) = [];

%% Commit the updated cell-vertex lists.
C{cu} = cell_list_cu;
C{cv} = cell_list_cv;
C{c1} = cell_list_c1;
C{c2} = cell_list_c2;

% --- Update vertex-vertex adjacency: break the old (u0-u2), (v0-v2)
% connections and form the new (u2-v0), (u0-v2) connections.
topo.VVadj(u0, u2) = 0; topo.VVadj(u2, u0) = 0;
topo.VVadj(v0, v2) = 0; topo.VVadj(v2, v0) = 0;
topo.VVadj(u2, v0) = 1; topo.VVadj(v0, u2) = 1;
topo.VVadj(u0, v2) = 1; topo.VVadj(v2, u0) = 1;

% --- Update cell-vertex adjacency to match the new cell vertex lists.
topo.CVadj(c1, u0) = 0;
topo.CVadj(c2, v0) = 0;
topo.CVadj(cu, v0) = 1;
topo.CVadj(cv, u0) = 1;

% --- Locate the edges that were rewired above (u0-u2 and v0-v2), so the
% edge list and vertex-edge / cell-edge adjacency can be updated too.
E1 = find(...
    (topo.E_list(:,1) == u0 & topo.E_list(:,2) == u2) | ...
    (topo.E_list(:,2) == u0 & topo.E_list(:,1) == u2));
E4 = find(...
    (topo.E_list(:,1) == v0 & topo.E_list(:,2) == v2) | ...
    (topo.E_list(:,2) == v0 & topo.E_list(:,1) == v2));

% Rewire those edges to their new endpoints.
topo.E_list(E1,:) = [v0, u2];
topo.E_list(E4,:) = [u0, v2];

% Update vertex-edge incidence for the rewired edges.
topo.VEadj(u0, E1) = 0; topo.VEadj(v0, E1) = 1;
topo.VEadj(v0, E4) = 0; topo.VEadj(u0, E4) = 1;

% Update cell-edge incidence: the flipped edge (T1_edge_ID) now belongs
% to cu/cv instead of c1/c2.
topo.CEadj(c1, T1_edge_ID) = 0;
topo.CEadj(c2, T1_edge_ID) = 0;
topo.CEadj(cu, T1_edge_ID) = 1;
topo.CEadj(cv, T1_edge_ID) = 1;

% --- Rotate the collapsed edge 90 degrees about its midpoint, so it
% re-expands perpendicular to its original orientation (this is what
% actually realizes the neighbor exchange geometrically). The edge
% vector is first unwrapped using the sheared minimum-image convention.
ex = V(v0,1) - V(u0,1);
ey = V(v0,2) - V(u0,2);
[ex, ey] = pbc_dist_LE(ex, ey, box_size, strain);

xv0 = V(u0,1) + ex;
yv0 = V(u0,2) + ey;
xu0 = V(u0,1);
yu0 = V(u0,2);
xmid = 0.5 * (xv0 + xu0);
ymid = 0.5 * (yv0 + yu0);

V(v0,1) = -(yv0-ymid) + xmid;
V(v0,2) =  (xv0-xmid) + ymid;
V(u0,1) = -(yu0-ymid) + xmid;
V(u0,2) =  (xu0-xmid) + ymid;

% Wrap the rotated vertices back into the box under sheared PBC.
V = position_mod_LE(V, box_size, strain);
end
