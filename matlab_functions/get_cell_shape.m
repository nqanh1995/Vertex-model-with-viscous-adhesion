% Periodic Boundary Conditions version
function [perim_list, area_list] = get_cell_shape(V, newC)
% GET_CELL_SHAPE  Compute perimeter and area of every cell polygon under
% periodic boundary conditions (PBC).
%
%   [perim_list, area_list] = GET_CELL_SHAPE(V, newC) computes, for each
%   cell, the perimeter and area of the polygon formed by its vertices,
%   correctly "unwrapping" vertices across the periodic box boundary
%   before applying the shoelace formula.
%
% Inputs:
%   V     - Nvert x 2 array of vertex positions
%   newC  - Ncell x 1 cell array; newC{i} lists the vertex indices (in
%           order around the polygon) belonging to cell i
%
% Outputs:
%   perim_list - Ncell x 1 vector of cell perimeters
%   area_list  - Ncell x 1 vector of cell areas

Npts = size(newC, 1);
box_size = sqrt(Npts);              % assumes a square box of unit cell density
area_list = zeros(Npts, 1);
perim_list = zeros(Npts, 1);

for i = 1:Npts
    X = V(newC{i}, :);   % vertex coordinates of cell i, in polygon order
    Z = size(X, 1);

    % Unwrap all vertices relative to the first vertex, so the polygon is
    % not artificially split by the periodic boundary.
    Ximg = zeros(Z, 2);
    Ximg(1,:) = X(1,:);
    Ximg(2:Z,1) = X(1,1) + pbc_dist(X(2:Z,1) - X(1,1), box_size);
    Ximg(2:Z,2) = X(1,2) + pbc_dist(X(2:Z,2) - X(1,2), box_size);

    % Edge vectors between consecutive (unwrapped) vertices.
    dx = Ximg([2:Z 1], 1) - Ximg(:,1);
    dy = Ximg([2:Z 1], 2) - Ximg(:,2);

    % Shoelace formula for polygon area (summation form of the boundary
    % integral); abs() makes the result orientation-independent.
    area_list(i) = abs(sum(Ximg(:,2).*dx - Ximg(:,1).*dy) / 2);

    % Perimeter = sum of edge lengths.
    perim_list(i) = sum(sqrt(dx.*dx + dy.*dy));
end
end
