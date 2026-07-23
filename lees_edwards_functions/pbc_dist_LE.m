function [dx, dy] = pbc_dist_LE(dx_old, dy_old, box_L, strain_value)
% PBC_DIST_LE  Minimum-image displacement under Lees-Edwards (sheared)
% periodic boundary conditions.
%
%   [dx, dy] = PBC_DIST_LE(dx_old, dy_old, box_L, strain_value) applies
%   the Lees-Edwards convention: when a displacement wraps across the top
%   or bottom of the box, it is additionally shifted horizontally by the
%   accumulated shear strain, before the usual minimum-image wrapping is
%   applied in x.
%
% Reference: http://homepage.univie.ac.at/franz.vesely/simsp/dx/node40.html
%
% Inputs:
%   dx_old, dy_old - raw displacement(s) before PBC correction
%   box_L          - size of the (square) periodic simulation box
%   strain_value   - current Lees-Edwards shear strain (image shift per
%                    box height, in units of box_L)
%
% Outputs:
%   dx, dy - displacement(s) wrapped to the nearest sheared periodic image

% Number of box heights the y-displacement wraps across; this determines
% how many "sheared images" the point has crossed.
cy = round(dy_old / box_L);

% Shift x by the shear accumulated over those wraps, then apply the
% standard minimum-image wrapping in x.
cdelx = dx_old - cy * strain_value * box_L;
dx = cdelx - round(cdelx / box_L) * box_L;

% Standard minimum-image wrapping in y.
dy = dy_old - cy * box_L;
end
