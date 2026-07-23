function ff = draw_state_LE(V, C, strain, fignum, c, lw, opac)
% DRAW_STATE_LE  Plot the tissue's cell polygons under Lees-Edwards
% sheared periodic boundary conditions.
%
%   ff = DRAW_STATE_LE(V, C, strain, fignum, c, lw, opac) draws the
%   outline of every cell, unwrapping each cell's vertices (relative to
%   its first vertex) using the sheared minimum-image convention so cells
%   that straddle the periodic boundary are drawn as a single contiguous
%   polygon.
%
% Inputs:
%   V      - Nvert x 2 array of vertex positions
%   C      - Ncell x 1 cell array; C{i} lists the vertex indices
%            belonging to cell i
%   strain - current Lees-Edwards shear strain
%   fignum - figure number to draw into (use 0 to create a new figure)
%   c      - (optional) line color, default 'k'
%   lw     - (optional) line width, default 1
%   opac   - (optional) line opacity/alpha, default 1
%
% Output:
%   ff - handle to the figure used
%
% Note: c, lw, opac are given defaults below (they were previously
% required positional arguments, causing an error whenever this function
% was called with just 4 arguments, as it is elsewhere in this codebase).
if nargin < 5 || isempty(c),    c = 'k';   end
if nargin < 6 || isempty(lw),   lw = 1;    end
if nargin < 7 || isempty(opac), opac = 1;  end

if fignum == 0
    ff = figure();
else
    ff = figure(fignum);
end

N = size(C, 1);
box_size = sqrt(N);

cm = colormap(jet(N)); %#ok<NASGU> % reserved for optional per-cell coloring
hold on;
for i = 1:N
    X = V(C{i}, 1);
    Y = V(C{i}, 2);
    z = size(X, 1);

    % Unwrap this cell's vertices relative to its first vertex, so the
    % polygon isn't artificially split by the sheared periodic boundary.
    dx = X(2:z) - X(1);
    dy = Y(2:z) - Y(1);
    [dx, dy] = pbc_dist_LE(dx, dy, box_size, strain);
    X(2:z) = X(1) + dx;
    Y(2:z) = Y(1) + dy;

    % Draw the closed polygon outline (repeat first vertex at the end).
    plot([X(:); X(1)], [Y(:); Y(1)], color=c, LineWidth=lw);
    alpha(opac);

    hold on;
end

axis([-2, box_size+2, -2, box_size+2])
rectangle('Position', [0, 0, box_size, box_size])   % draw the simulation box
pbaspect([1, 1, 1]);
set(gca, 'Color', 'none');   % transparent axes background
end
