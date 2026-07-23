function [new_x, new_y] = drifting_correction(x, y)
% DRIFTING_CORRECTION  Remove uniform (whole-system) translational drift
% from a time series of positions.
%
%   [new_x, new_y] = DRIFTING_CORRECTION(x, y) takes position trajectories
%   x, y (rows = time steps, columns = individual points/vertices) and
%   subtracts, at every time step, the mean displacement averaged over all
%   points. This removes any net drift of the whole system (e.g. due to
%   numerical or boundary artifacts) while preserving relative motion
%   between points.
%
% Inputs:
%   x, y  - position trajectories, size (Nsteps x Npoints)
%
% Output:
%   new_x, new_y - drift-corrected trajectories, same size as x, y

% Step-to-step displacements for every point.
global_dx = diff(x, 1, 1);
global_dy = diff(y, 1, 1);

% Average displacement across all points at each step = the drift.
mean_global_dx = mean(global_dx, 2);
mean_global_dy = mean(global_dy, 2);

% Remove the drift from each point's displacement.
non_drift_dx = bsxfun(@minus, global_dx, mean_global_dx);
non_drift_dy = bsxfun(@minus, global_dy, mean_global_dy);

% Reconstruct drift-free trajectories by integrating the corrected
% displacements back up, starting from the original initial positions.
new_x = bsxfun(@plus, cumsum(non_drift_dx), x(1,:));
new_y = bsxfun(@plus, cumsum(non_drift_dy), y(1,:));
end
