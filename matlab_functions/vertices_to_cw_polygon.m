% order a point set such that the points are in counter clockwise order
% from the centeroid
function[order] = vertices_to_cw_polygon(x,y);
cx = mean(x);
cy = mean(y);
a = atan2(y - cy, x - cx);
[~, order] = sort(a);
% newx = x(order);
% newy = y(order);