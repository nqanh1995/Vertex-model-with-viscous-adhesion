function [stress_tensor, cell_stress] = get_viscous_qvm_stress(V, C, topo, KP, P0, KA, A0, strain, Fvx, Fvy)
% GET_VISCOUS_QVM_STRESS  Compute the tissue-level and per-cell stress
% tensors, including the viscous (cell-cell adhesion) force contribution,
% under Lees-Edwards (sheared) periodic boundary conditions.
%
%   [stress_tensor, cell_stress] = GET_VISCOUS_QVM_STRESS(V, C, topo, KP,
%   P0, KA, A0, strain, Fvx, Fvy) computes each cell's contribution to
%   the stress tensor from three sources: (1) isotropic pressure from
%   area elasticity, (2) line tension from perimeter elasticity along
%   each edge, and (3) the symmetrized virial of the viscous coupling
%   force (Fvx, Fvy) acting at each vertex relative to its cell's
%   centroid. Cell-level stresses are then area-weighted and summed to
%   give the intensive tissue-level stress tensor.
%
% Inputs:
%   V, C     - vertex positions and cell-vertex lists
%   topo     - topology struct (edge list, cell-edge adjacency, counts)
%   KP, P0   - perimeter elasticity coefficient and target perimeter, per
%              cell
%   KA, A0   - area elasticity coefficient and target area, per cell
%   strain   - current Lees-Edwards shear strain
%   Fvx, Fvy - Nvert x 1 vectors of the viscous coupling force at each
%              vertex (as computed in get_vertex_velo_viscous.m)
%
% Outputs:
%   stress_tensor - 2x2 tissue-level (intensive) stress tensor
%   cell_stress   - Ncell x 4 per-cell stress components, ordered as
%                   [sigma_xx, sigma_xy, sigma_yx, sigma_yy]

[cell_perimeters, cell_areas] = get_cell_shape_LE(V, C, strain);

% Line tension on every edge: 2*KP*(perimeter - P0) for the cell, spread
% across all of that cell's edges via the cell-edge adjacency matrix.
tensions = (2*KP.*(cell_perimeters - P0))' * topo.CEadj;
tensions = tensions(:);

% Isotropic pressure in every cell: -2*KA*(area - A0).
cell_pressures = -2*KA.*(cell_areas-A0);
cell_pressures = cell_pressures(:);

% Per-cell stress tensor components, ordered as sxx, sxy, syx, syy.
cell_stress = zeros(topo.Ncell, 4);

box_size = sqrt(topo.Ncell);

% Edge vectors and lengths, corrected for Lees-Edwards sheared PBC.
ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
[ex, ey] = pbc_dist_LE(ex, ey, box_size, strain);
edge_L = sqrt(ex.^2 + ey.^2);

stress_tensor = zeros(2, 2);

% --- Pressure contribution: isotropic, on the diagonal only.
for c = 1:topo.Ncell
    cell_stress(c,1) = -cell_pressures(c);
    cell_stress(c,4) = -cell_pressures(c);
end

% --- Tension contribution: each edge of a cell adds a
% tension * (edge-vector outer product) / edge-length term, normalized by
% twice the cell area.
for c = 1:topo.Ncell
    eid = topo.CEadj(c,:);
    cell_stress(c,1) = cell_stress(c,1) + ...
        sum(tensions(eid).*ex(eid).*ex(eid)./edge_L(eid)) / (2*cell_areas(c));
    cell_stress(c,2) = cell_stress(c,2) + ...
        sum(tensions(eid).*ex(eid).*ey(eid)./edge_L(eid)) / (2*cell_areas(c));
    cell_stress(c,3) = cell_stress(c,3) + ...
        sum(tensions(eid).*ey(eid).*ex(eid)./edge_L(eid)) / (2*cell_areas(c));
    cell_stress(c,4) = cell_stress(c,4) + ...
        sum(tensions(eid).*ey(eid).*ey(eid)./edge_L(eid)) / (2*cell_areas(c));
end

% --- Viscous force contribution: symmetrized virial stress
%   sigma_ab += -sum(R_a*F_b + F_a*R_b) / (2*z*Area)
% where R is each vertex's position relative to the cell centroid and F
% is the viscous coupling force on that vertex (Fvx, Fvy).
for c = 1:topo.Ncell
    verts = C{c};
    z = length(verts);
    Vx = V(verts,1);
    Vy = V(verts,2);

    % Unwrap this cell's vertices relative to its first vertex (sheared
    % PBC), then compute the centroid and each vertex's offset from it.
    [dvx, dvy] = pbc_dist_LE(Vx(2:z)-Vx(1), Vy(2:z)-Vy(1), box_size, strain);
    Vx(2:z) = Vx(1) + dvx;
    Vy(2:z) = Vy(1) + dvy;
    ccx = mean(Vx, 1);
    ccy = mean(Vy, 1);
    [Rx, Ry] = pbc_dist_LE(Vx-ccx, Vy-ccy, box_size, strain);

    cell_stress(c,1) = cell_stress(c,1) - sum(Rx.*Fvx(verts) + Rx.*Fvx(verts)) / (2*z*cell_areas(c));
    cell_stress(c,2) = cell_stress(c,2) - sum(Rx.*Fvy(verts) + Ry.*Fvx(verts)) / (2*z*cell_areas(c));
    cell_stress(c,3) = cell_stress(c,3) - sum(Ry.*Fvx(verts) + Rx.*Fvy(verts)) / (2*z*cell_areas(c));
    cell_stress(c,4) = cell_stress(c,4) - sum(Ry.*Fvy(verts) + Ry.*Fvy(verts)) / (2*z*cell_areas(c));
end

% --- Combine cell-level stresses into the tissue-level stress tensor
% (area-weighted sum), then normalize by total area to make it intensive.
stress_tensor(1,1) = sum(cell_stress(:,1).*cell_areas);
stress_tensor(1,2) = sum(cell_stress(:,2).*cell_areas);
stress_tensor(2,1) = sum(cell_stress(:,3).*cell_areas);
stress_tensor(2,2) = sum(cell_stress(:,4).*cell_areas);

stress_tensor = stress_tensor / sum(cell_areas);
end
