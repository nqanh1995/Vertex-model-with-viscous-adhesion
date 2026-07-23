function V = position_mod_LE(V, box_size, strain)
% POSITION_MOD_LE  Wrap vertex positions back into the primary simulation
% box under Lees-Edwards (sheared) periodic boundary conditions.
%
%   V = POSITION_MOD_LE(V, box_size, strain) wraps y-coordinates into
%   [0, box_size) normally, and wraps x-coordinates into [0, box_size),
%   additionally applying the Lees-Edwards horizontal shear shift for any
%   vertex that crossed the top or bottom of the box (since, under shear,
%   the periodic image above/below the box is offset horizontally by
%   strain*box_size).
%
% Inputs:
%   V        - Nvert x 2 array of vertex positions
%   box_size - size of the (square) periodic simulation box
%   strain   - current Lees-Edwards shear strain
%
% Output:
%   V - vertex positions wrapped into the primary box

Vx = V(:,1);
Vy = V(:,2);

% Vertex crossed above the top of the box: shift x by the shear offset
% before wrapping (entering from the bottom, sheared image).
ind = Vy >= box_size;
Vx(ind) = mod(Vx(ind) - strain*box_size, box_size);

% Vertex already inside the box in y: just wrap x normally.
ind = Vy >= 0 & Vy < box_size;
Vx(ind) = mod(Vx(ind), box_size);

% Vertex crossed below the bottom of the box: shift x by the opposite
% shear offset before wrapping (entering from the top, sheared image).
ind = Vy < 0;
Vx(ind) = mod(Vx(ind) + strain*box_size, box_size);

% Wrap y normally.
Vy = mod(Vy, box_size);

V = [Vx, Vy];
end
