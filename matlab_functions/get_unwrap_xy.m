function [xcm_u, ycm_u] = get_unwrap_xy(xcm, ycm)
% GET_UNWRAP_XY  Unwrap center-of-mass trajectories across periodic
% boundary crossings, so trajectories are continuous instead of jumping
% by a box length whenever a particle crosses the periodic boundary.
%
%   [xcm_u, ycm_u] = GET_UNWRAP_XY(xcm, ycm) takes trajectory matrices
%   (rows = time steps, columns = particles) and removes the artificial
%   jumps caused by periodic wrapping, producing "unwrapped" trajectories
%   suitable for computing displacements, MSDs, etc.
%
% Inputs:
%   xcm, ycm - T x npts matrices of (possibly wrapped) positions
%
% Outputs:
%   xcm_u, ycm_u - T x npts matrices of unwrapped positions

[T, npts] = size(xcm);
box_size = sqrt(npts);   % assumes a square box of unit cell density

xcm_u = zeros(T, npts);
ycm_u = zeros(T, npts);

for j = 1:npts
    xcm_u(:,j) = unwrap_1d(xcm(:,j), box_size, T);
    ycm_u(:,j) = unwrap_1d(ycm(:,j), box_size, T);
end
end

function x_u = unwrap_1d(x, box_size, T)
% Unwrap a single 1D trajectory: detect steps larger than half the box
% (a periodic wrap) and accumulate a compensating offset from that point
% onward.
boundary_pts = (abs(diff(x)) > box_size/2) .* sign(diff(x));
boundary_idx = find(boundary_pts ~= 0);

wraplist = zeros(T, 1);
for i = 1:size(boundary_idx, 1)
    idx = boundary_idx(i);
    idx_sign = boundary_pts(idx);
    wraplist(idx+1:end) = wraplist(idx+1:end) - idx_sign * box_size;
end

x_u = x + wraplist;
end
