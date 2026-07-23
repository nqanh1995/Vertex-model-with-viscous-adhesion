% Periodic Boundary Conditions version
function [perimeters, areas, a1, a2, orientation] = cell_shape_info(V, newC)
% CELL_SHAPE_INFO  Compute shape descriptors (area, perimeter, principal
% axes, orientation) for every cell polygon under periodic boundary
% conditions (PBC).
%
%   [perimeters, areas, a1, a2, orientation] = CELL_SHAPE_INFO(V, newC)
%   unwraps each cell's vertices across the periodic box, orders them
%   around the polygon, and passes them to POLYGEOM to obtain geometric
%   and principal-moment-of-area properties.
%
% Inputs:
%   V     - Nvert x 2 array of vertex positions
%   newC  - Ncell x 1 cell array; newC{i} lists the vertex indices
%           belonging to cell i
%
% Outputs:
%   perimeters  - Ncell x 1 vector of cell perimeters
%   areas       - Ncell x 1 vector of cell areas
%   a1, a2      - Ncell x 1 vectors of principal semi-axis lengths
%                 (from the cell's principal moments of area)
%   orientation - Ncell x 1 vector of cell orientation angles
%
% Note: requires the third-party POLYGEOM function (MATLAB File
% Exchange) to be on the path; it is not included in this repository.

Npts = size(newC, 1);
box_size = sqrt(Npts);              % assumes a square box of unit cell density
areas = zeros(Npts, 1);
perimeters = zeros(Npts, 1);
orientation = zeros(Npts, 1);
a1 = zeros(Npts, 1);
a2 = zeros(Npts, 1);

for i = 1:Npts
    vx = V(newC{i}, 1);
    vy = V(newC{i}, 2);
    z = size(vx, 1);

    % Unwrap vertices relative to the first one, so the polygon is not
    % artificially split by the periodic boundary.
    vx(2:z) = vx(1) + pbc_dist(vx(2:z) - vx(1), box_size);
    vy(2:z) = vy(1) + pbc_dist(vy(2:z) - vy(1), box_size);

    % Re-order vertices sequentially around the polygon (required by
    % polygeom, which assumes a properly ordered boundary).
    ccw_order = vertices_to_cw_polygon(vx, vy);
    vx = vx(ccw_order);
    vy = vy(ccw_order);

    % geom: [area, centroid_x, centroid_y, perimeter]
    % cpmo: centroidal principal moments of area -> [Ixx, angle, Iyy] etc.
    [geom, ~, cpmo] = polygeom(vx, vy);
    areas(i) = geom(1);
    perimeters(i) = geom(4);
    a1(i) = cpmo(3);
    a2(i) = cpmo(1);
    orientation(i) = cpmo(2);
end
end
