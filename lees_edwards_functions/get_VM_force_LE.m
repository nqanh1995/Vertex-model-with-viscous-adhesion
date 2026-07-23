function [Vertex_Force, Cell_Perim, Cell_Area] = get_VM_force_LE(V, C, KP, P0, KA, A0, strain)
% GET_VM_FORCE_LE  Compute the net mechanical force on every vertex from
% the standard vertex-model energy (area + perimeter elasticity), under
% Lees-Edwards (sheared) periodic boundary conditions.
%
%   [Vertex_Force, Cell_Perim, Cell_Area] = GET_VM_FORCE_LE(V, C, KP, P0,
%   KA, A0, strain) computes, for every cell, the analytical gradient of
%   its energy
%       E_i = KA_i*(Area_i - A0_i)^2 + KP_i*(Perim_i - P0_i)^2
%   with respect to each of its vertices' positions, and assembles the
%   force on each vertex as (minus) the sum of these gradient
%   contributions from all cells sharing that vertex.
%
% Inputs:
%   V      - Nvert x 2 array of vertex positions
%   C      - Ncell x 1 cell array of vertex indices (in polygon order)
%            belonging to each cell
%   KP, P0 - perimeter elasticity coefficient and target perimeter, per
%            cell
%   KA, A0 - area elasticity coefficient and target area, per cell
%   strain - current Lees-Edwards shear strain
%
% Outputs:
%   Vertex_Force - Nvert x 2 array of net force on each vertex
%   Cell_Perim   - Ncell x 1 vector of current cell perimeters
%   Cell_Area    - Ncell x 1 vector of current cell areas

V_count = size(V, 1);
N = length(C);
box_size = sqrt(N);

[Cell_Perim, Cell_Area] = get_cell_shape_LE(V, C, strain);

% Accumulator for each cell's contribution to the force on every vertex:
% row i holds cell i's gradient contributions, spread across the columns
% for the vertices it touches (x-components in the first V_count
% columns, y-components in the last V_count columns).
dE_H_mat = zeros(N, 2*V_count);

for i = 1:N   % loop over cells
    neibs = C{i};
    z = length(neibs);
    Vx = V(neibs, 1);
    Vy = V(neibs, 2);

    % Unwrap this cell's vertices relative to its last vertex, so all
    % vertices are expressed in the same (sheared) periodic image before
    % computing finite differences between them.
    [dx, dy] = pbc_dist_LE(Vx - Vx(end), Vy - Vy(end), box_size, strain);
    Vx = Vx(end) + dx;
    Vy = Vy(end) + dy;

    % Extend the vertex list cyclically so each vertex has a well-defined
    % previous/next neighbor around the polygon, including wrap-around.
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];

    for n = 2:z+1
        % Gradient of the area term (shoelace formula) w.r.t. vertex n.
        Ehx = KA(i) * (Cell_Area(i)-A0(i)) * (Vy(n+1) - Vy(n-1));
        Ehy = KA(i) * (Cell_Area(i)-A0(i)) * (Vx(n-1) - Vx(n+1));

        % Gradient of the perimeter term w.r.t. vertex n, from the two
        % edges meeting at that vertex (to the previous and next vertex).
        minus_x = 2*KP(i)*(Cell_Perim(i)-P0(i)) * ...
            (Vx(n)-Vx(n-1)) / sqrt((Vx(n-1)-Vx(n))^2 + (Vy(n-1)-Vy(n))^2);
        plus_x = 2*KP(i)*(Cell_Perim(i)-P0(i)) * ...
            (Vx(n)-Vx(n+1)) / sqrt((Vx(n+1)-Vx(n))^2 + (Vy(n+1)-Vy(n))^2);

        minus_y = 2*KP(i)*(Cell_Perim(i)-P0(i)) * ...
            (Vy(n)-Vy(n-1)) / sqrt((Vx(n-1)-Vx(n))^2 + (Vy(n-1)-Vy(n))^2);
        plus_y = 2*KP(i)*(Cell_Perim(i)-P0(i)) * ...
            (Vy(n)-Vy(n+1)) / sqrt((Vx(n+1)-Vx(n))^2 + (Vy(n+1)-Vy(n))^2);

        % Guard against zero-length edges (division by zero -> NaN):
        % simply omit that edge's contribution rather than propagating NaN.
        if ~isnan(minus_x) && ~isnan(minus_y)
            Ehx = Ehx + minus_x;
            Ehy = Ehy + minus_y;
        end
        if ~isnan(plus_x) && ~isnan(plus_y)
            Ehx = Ehx + plus_x;
            Ehy = Ehy + plus_y;
        end

        dE_H_mat(i, neibs(n)) = Ehx;
        dE_H_mat(i, neibs(n)+V_count) = Ehy;
    end
end

% Force = -dE/dV, summed over every cell's contribution to each vertex.
Vertex_Force = -sum(dE_H_mat, 1);
Vertex_Force = [Vertex_Force(1:V_count); Vertex_Force(V_count+1:end)]';
end
