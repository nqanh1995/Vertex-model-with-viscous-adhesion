function [vvx, vvy] = get_vertex_velo_viscous(V, topo, VF, B, strain)
% GET_VERTEX_VELO_VISCOUS  Solve for vertex velocities including
% cell-cell viscous adhesion coupling between neighboring vertices
% (Lees-Edwards sheared PBC version).
%
%   [vvx, vvy] = GET_VERTEX_VELO_VISCOUS(V, topo, VF, B, strain) builds
%   the linear system that couples each vertex's velocity to its
%   neighbors' velocities through a viscous drag term (coefficient B),
%   in addition to the usual local (mechanical + active) force VF, and
%   solves for the resulting vertex velocities.
%
% Inputs:
%   V      - Nvert x 2 array of vertex positions
%   topo   - topology struct (must contain VVadj, the vertex-vertex
%            adjacency matrix)
%   VF     - Nvert x 2 array of net force on each vertex (mechanical +
%            active), before viscous coupling
%   B      - cell-cell viscous coefficient (xi)
%   strain - current Lees-Edwards shear strain
%
% Outputs:
%   vvx, vvy - Nvert x 1 vectors of vertex velocities

N = size(V, 1) / 2;
box_size = sqrt(N);
adjVV = topo.VVadj;

Vx = V(:,1);
Vy = V(:,2);

% Neighbor coordinates broadcast through the adjacency mask, used to
% build the pairwise displacement (Vi - Vj) matrices between adjacent
% vertices.
neib_Vx = adjVV .* Vx';
neib_Vy = adjVV .* Vy';
xi_xj = -(neib_Vx - Vx .* adjVV);   % Vxi - Vxj matrix, j neighboring i
yi_yj = -(neib_Vy - Vy .* adjVV);   % Vyi - Vyj matrix, j neighboring i

% Apply the minimum-image convention under Lees-Edwards shear, then flip
% sign back to get the (i - j) displacement convention.
[ex, ey] = pbc_dist_LE(xi_xj, yi_yj, box_size, strain);
ex = -ex;
ey = -ey;

% Edge length matrix (matches the adjVV sparsity pattern) and unit
% direction vectors along each vertex-vertex edge.
edge_L = sqrt(ex.^2 + ey.^2);
ex_hat = ex ./ edge_L;
ex_hat(isnan(ex_hat)) = 0;   % non-neighbor entries (0/0) -> 0
ey_hat = ey ./ edge_L;
ey_hat(isnan(ey_hat)) = 0;

% Build the coefficient blocks of the viscous coupling system: diagonal
% terms are "self drag" (identity + neighbor viscous contributions),
% off-diagonal terms are "coupling drag" to each neighbor.
coef_xx = eye(size(adjVV)) .* (B * sum(ex.*ex_hat, 2) + 1) - B * ex .* ex_hat;
coef_xy = eye(size(adjVV)) .* (B * sum(ey.*ex_hat, 2))     - B * ey .* ex_hat;
coef_yx = eye(size(adjVV)) .* (B * sum(ey.*ex_hat, 2))     - B * ex .* ey_hat;
coef_yy = eye(size(adjVV)) .* (B * sum(ey.*ey_hat, 2) + 1) - B * ey .* ey_hat;

% Assemble the full block system [x-eqns; y-eqns] coupling all vertices'
% x and y velocity components together.
coeff_mat = zeros(4*N, 4*N);
coeff_mat(1:2*N, 1:2*N)         = coef_xx;
coeff_mat(1:2*N, 2*N+1:4*N)     = coef_xy;
coeff_mat(2*N+1:4*N, 1:2*N)     = coef_yx;
coeff_mat(2*N+1:4*N, 2*N+1:4*N) = coef_yy;

A = [VF(:,1); VF(:,2)];

% Solve the linear system for vertex velocities (backslash is faster and
% more numerically stable than forming the explicit matrix inverse).
vv = coeff_mat \ A;
vvx = vv(1:2*N);
vvy = vv(2*N+1:end);
end
