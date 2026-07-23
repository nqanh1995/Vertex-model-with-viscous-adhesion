function ff = draw_edges_LE(V, topo, strain, fignum)
% DRAW_EDGES_LE  Plot every edge of the tissue network (rather than cell
% polygons), each drawn in a distinct color, under Lees-Edwards sheared
% PBC.
%
%   ff = DRAW_EDGES_LE(V, topo, strain, fignum) draws each edge as a
%   short line segment from its first vertex, using the sheared
%   minimum-image displacement to its second vertex (so edges crossing
%   the periodic boundary are drawn correctly rather than stretching
%   across the whole box). Useful for visually inspecting topology/T1
%   transitions.
%
% Inputs:
%   V      - Nvert x 2 array of vertex positions
%   topo   - topology struct (must contain Ncell, Nedge, E_list)
%   strain - current Lees-Edwards shear strain
%   fignum - figure number to draw into (use 0 to create a new figure)
%
% Output:
%   ff - handle to the figure used

if fignum == 0
    ff = figure();
else
    ff = figure(fignum);
end

box_size = sqrt(topo.Ncell);

% Give each edge a distinct color so individual edges are easy to track
% across frames (e.g. to watch a specific T1 transition).
cm = colormap(hsv(topo.Nedge));
hold on;
for i = 1:topo.Nedge
    v1 = topo.E_list(i,1);
    v2 = topo.E_list(i,2);
    dx = V(v2,1) - V(v1,1);
    dy = V(v2,2) - V(v1,2);
    [dx, dy] = pbc_dist_LE(dx, dy, box_size, strain);

    plot([V(v1,1), V(v1,1)+dx], [V(v1,2), V(v1,2)+dy], '-', 'color', cm(i,:), 'LineWidth', 2);
end

axis([-2, box_size+2, -2, box_size+2])
rectangle('Position', [0, 0, box_size, box_size])   % draw the simulation box
pbaspect([1, 1, 1]);
set(gca, 'Color', 'none');   % transparent axes background
text(box_size/2, -1, ['strain = ', num2str(strain)])
end
