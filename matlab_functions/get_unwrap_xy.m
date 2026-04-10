function [xcm_u, ycm_u] = get_unwrap_xy(xcm,ycm)
[T, npts] = size(xcm);
box_size = sqrt(npts);

% Unwrapping x and y positions
xcm_u = zeros(T,npts);
ycm_u = zeros(T,npts);

for j = 1:npts
    boundary_pts=(abs(diff(xcm(:,j)))>box_size/2) .* sign(diff(xcm(:,j)));
    boundary_idx = find(boundary_pts~=0);
    wraplist = zeros(T,1);
    for i = 1:size(boundary_idx,1)
        idx = boundary_idx(i);
        idx_sign = boundary_pts(idx);
        wraplist(idx+1:end) = wraplist(idx+1:end)-idx_sign*box_size;

    end
    xcm_u(:,j) = xcm(:,j)+wraplist;

    boundary_pts=(abs(diff(ycm(:,j)))>box_size/2) .* sign(diff(ycm(:,j)));
    boundary_idx = find(boundary_pts~=0);
    wraplist = zeros(T,1);
    for i = 1:size(boundary_idx,1)
        idx = boundary_idx(i);
        idx_sign = boundary_pts(idx);
        wraplist(idx+1:end) = wraplist(idx+1:end)-idx_sign*box_size;

    end
    ycm_u(:,j) = ycm(:,j)+wraplist;
end
end