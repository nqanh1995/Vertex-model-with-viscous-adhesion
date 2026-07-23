function ff = draw_state(V, C)
% DRAW_STATE  Plot the tissue's cell polygons (non-sheared PBC version).
%
%   ff = DRAW_STATE(V, C) draws the outline of every cell, unwrapped
%   across the periodic boundary so cells that straddle the box edge are
%   drawn as a single contiguous polygon rather than being split.
%
% Inputs:
%   V - Nvert x 2 array of vertex positions
%   C - Ncell x 1 cell array; C{i} lists the vertex indices belonging to
%       cell i
%
% Output:
%   ff - handle to the created figure

ff = figure();
N = size(C, 1);
box_size = sqrt(N);

hold on;
for i = 1:N
    X = V(C{i}, :);
    Z = size(X, 1);

    % Unwrap this cell's vertices relative to its first vertex, so the
    % polygon isn't artificially split by the periodic boundary.
    Ximg = zeros(Z, 2);
    Ximg(1,:) = X(1,:);
    Ximg(2:Z,1) = X(1,1) + pbc_dist(X(2:Z,1) - X(1,1), box_size);
    Ximg(2:Z,2) = X(1,2) + pbc_dist(X(2:Z,2) - X(1,2), box_size);

    % If the unwrapped polygon's centroid ends up on the "wrong side" of
    % the box relative to the reference vertex, shift the whole polygon
    % back by one box length so it's drawn inside/near the visible box.
    meanxy = mean(Ximg, 1);
    if meanxy(1) - X(1,1) > 0.5*box_size
        Ximg(:,1) = Ximg(:,1) - box_size;
    end
    if meanxy(1) - X(1,1) < -0.5*box_size
        Ximg(:,1) = Ximg(:,1) + box_size;
    end
    if meanxy(2) - X(1,2) > 0.5*box_size
        Ximg(:,2) = Ximg(:,2) - box_size;
    end
    if meanxy(2) - X(1,2) < -0.5*box_size
        Ximg(:,2) = Ximg(:,2) + box_size;
    end

    % Draw the closed polygon outline (repeat first vertex at the end).
    plot([Ximg(:,1); Ximg(1,1)], [Ximg(:,2); Ximg(1,2)], '-k');
    hold on;
end

axis([-2, box_size+2, -2, box_size+2])
rectangle('Position', [0, 0, box_size, box_size])   % draw the simulation box
pbaspect([1, 1, 1]);
set(gca, 'Color', 'none');   % transparent axes background
end
