function new_dist = pbc_dist(old_dist, box_size)
% PBC_DIST  Apply the minimum-image convention under periodic boundary
% conditions (PBC) to a displacement (or array of displacements).
%
%   new_dist = PBC_DIST(old_dist, box_size) wraps any displacement whose
%   magnitude exceeds half the box size back into the range
%   [-box_size/2, box_size/2], which is the standard "minimum image"
%   correction used when computing distances in a periodic box.
%
% Inputs:
%   old_dist  - raw displacement(s) (scalar or array)
%   box_size  - size of the (square) periodic simulation box
%
% Output:
%   new_dist  - displacement(s) wrapped to the nearest periodic image

% Flag entries whose displacement is more than half a box width: these
% are "too far" and actually correspond to the neighboring periodic image.
pbc_switch = abs(old_dist) > 0.5 * box_size;

% For flagged entries, subtract one box length in the direction of the
% displacement (bringing it to the nearest image); unflagged entries are
% left unchanged.
new_dist = (1 - pbc_switch) .* old_dist + pbc_switch .* (old_dist - sign(old_dist) .* box_size);
return
