function [theta_v, edge_L] = get_theta_vertex(topo, V, theta, strain)
% GET_THETA_VERTEX  Map cell-level polarity angles onto vertex-level unit
% polarity vectors, using an edge-length-weighted average of the
% surrounding cells' polarities (Lees-Edwards sheared PBC version).
%
%   [theta_v, edge_L] = GET_THETA_VERTEX(topo, V, theta, strain) computes,
%   for every vertex, a weighted average of the polarity directions of
%   the cells meeting at that vertex, weighted by the length of the
%   shared edges, then normalizes the result to a unit vector. This is
%   used to project each cell's active self-propulsion direction onto its
%   vertices for the equations of motion.
%
% Inputs:
%   topo   - topology struct (adjacency matrices, edge list, counts)
%   V      - Nvert x 2 array of vertex positions
%   theta  - Ncell x 1 vector of cell polarity angles
%   strain - current Lees-Edwards shear strain
%
% Outputs:
%   theta_v - Nvert x 2 array of unit polarity vectors at each vertex
%   edge_L  - Nedge x 1 vector of edge lengths (sheared PBC minimum image)

box_size = sqrt(topo.Ncell);
Ncell = topo.Ncell; %#ok<NASGU> % kept for readability / potential future use
Nvert = topo.Nvert;

% Edge vectors and lengths, corrected for Lees-Edwards sheared PBC.
edge_dx = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
edge_dy = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
[edge_dx, edge_dy] = pbc_dist_LE(edge_dx, edge_dy, box_size, strain);
edge_L = sqrt(edge_dx.^2 + edge_dy.^2);

% Total perimeter contribution of each vertex (sum of adjacent edge
% lengths), used to normalize the weighting below.
L = topo.VEadj * edge_L;

% Unit polarity vector for every cell.
theta_c = [cos(theta), sin(theta)];

theta_v = zeros(Nvert, 2);
for i = 1:Nvert
    % Cells that touch vertex i.
    neibs_C = find(topo.CVadj(:,i));
    z = length(neibs_C);
    n_c = zeros(z, 2);
    n_c(:,1) = theta_c(neibs_C,1);
    n_c(:,2) = theta_c(neibs_C,2);
    a_c = zeros(z, 1);

    % Edges incident to vertex i.
    E_Vi = find(topo.VEadj(i,:));

    for j = 1:z
        % For each neighboring cell, find the edge(s) it shares with
        % vertex i, and weight its polarity by that shared edge length.
        cid = neibs_C(j);
        E_Cj = find(topo.CEadj(cid,:));
        l12 = intersect(E_Cj, E_Vi);
        a_c(j) = sum(edge_L(l12)) ./ (2 * z * L(i));
    end

    % Weighted sum of neighboring cell polarities.
    theta_v(i,:) = sum(a_c .* n_c, 1);

    % Normalize to a unit vector (direction only).
    magnitude = norm(theta_v(i,:));
    if magnitude > 0
        theta_v(i,:) = theta_v(i,:) / magnitude;
    else
        warning('theta_v at vertex %d has zero magnitude and cannot be normalized.', i);
    end
end
end
