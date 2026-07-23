% Periodic Boundary Conditions version
function [plist, alist] = get_cell_shape_LE(V, C, strain)
% GET_CELL_SHAPE_LE  Compute perimeter and area of every cell polygon
% under Lees-Edwards (sheared) periodic boundary conditions.
%
%   [plist, alist] = GET_CELL_SHAPE_LE(V, C, strain) computes, for each
%   cell, the perimeter and area of the polygon formed by its vertices,
%   unwrapping vertices across the sheared periodic boundary before
%   applying the shoelace formula.
%
% Inputs:
%   V      - Nvert x 2 array of vertex positions
%   C      - Ncell x 1 cell array; C{i} lists the vertex indices (in
%            order around the polygon) belonging to cell i
%   strain - current Lees-Edwards shear strain
%
% Outputs:
%   plist - Ncell x 1 vector of cell perimeters
%   alist - Ncell x 1 vector of cell areas

N = size(C, 1);
box_size = sqrt(N);
alist = zeros(N, 1);
plist = zeros(N, 1);

for i = 1:N
    X = V(C{i}, 1);
    Y = V(C{i}, 2);
    z = size(X, 1);

    % Unwrap all vertices relative to the first vertex, accounting for
    % Lees-Edwards shear, so the polygon is not artificially split by
    % the periodic boundary.
    dx = X(2:z) - X(1);
    dy = Y(2:z) - Y(1);
    [dx, dy] = pbc_dist_LE(dx, dy, box_size, strain);
    X(2:z) = X(1) + dx;
    Y(2:z) = Y(1) + dy;

    % Edge vectors between consecutive (unwrapped) vertices.
    dx = X([2:z, 1]) - X;
    dy = Y([2:z, 1]) - Y;

    % Shoelace formula for polygon area; abs() makes the result
    % orientation-independent.
    alist(i) = abs(sum(Y.*dx - X.*dy) / 2);

    % Perimeter = sum of edge lengths.
    plist(i) = sum(sqrt(dx.*dx + dy.*dy));
end
end
