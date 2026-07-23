function stress_tensor = calc_stress(V, C, topo, tensions, cell_pressures, strain)
% CALC_STRESS  Compute the tissue-level (Batchelor/virial) stress tensor
% from edge tensions and cell pressures, under Lees-Edwards (sheared)
% periodic boundary conditions.
%
%   stress_tensor = CALC_STRESS(V, C, topo, tensions, cell_pressures,
%   strain) combines an isotropic pressure contribution (from each cell's
%   area and pressure) with a line-tension contribution (from each edge's
%   tension and orientation) into the standard virial stress tensor for a
%   confluent tissue, then normalizes by total area to make it intensive.
%
% Inputs:
%   V, C           - vertex positions and cell-vertex lists
%   topo           - topology struct (edge list, edge count)
%   tensions       - Nedge x 1 vector of line tensions on each edge
%   cell_pressures - Ncell x 1 vector of pressure in each cell
%   strain         - current Lees-Edwards shear strain
%
% Output:
%   stress_tensor - 2x2 tissue-level stress tensor

[~, cell_areas] = get_cell_shape_LE(V, C, strain);

box_size = sqrt(topo.Ncell);

% Edge vectors and lengths, corrected for Lees-Edwards sheared PBC.
ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
[ex, ey] = pbc_dist_LE(ex, ey, box_size, strain);
edge_L = sqrt(ex.^2 + ey.^2);

stress_tensor = zeros(2, 2);

% Isotropic pressure contribution: -sum(pressure * area) on the diagonal.
stress_tensor(1,1) = -sum(cell_pressures .* cell_areas);
stress_tensor(2,2) = -sum(cell_pressures .* cell_areas);

% Line-tension contribution: each edge adds an outer-product-like term
% (tension * edge-vector-outer-product / edge-length) to the tensor.
for i = 1:topo.Nedge
    stress_tensor(1,1) = stress_tensor(1,1) + tensions(i)*ex(i)*ex(i)/edge_L(i);
    stress_tensor(2,2) = stress_tensor(2,2) + tensions(i)*ey(i)*ey(i)/edge_L(i);
    stress_tensor(1,2) = stress_tensor(1,2) + tensions(i)*ex(i)*ey(i)/edge_L(i);
    stress_tensor(2,1) = stress_tensor(2,1) + tensions(i)*ey(i)*ex(i)/edge_L(i);
end

% Normalize by total tissue area to make the stress an intensive quantity.
stress_tensor = stress_tensor / sum(cell_areas);
end
