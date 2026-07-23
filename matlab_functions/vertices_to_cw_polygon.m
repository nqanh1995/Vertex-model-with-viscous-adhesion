function [order] = vertices_to_cw_polygon(x, y)
% VERTICES_TO_CW_POLYGON  Sort a set of 2D points by angle around their
% centroid, so that they form a simple (non-self-intersecting) polygon.
%
%   order = VERTICES_TO_CW_POLYGON(x, y) returns the permutation "order"
%   such that x(order), y(order) traces the points in increasing angle
%   from the centroid (i.e. sequentially around the polygon, rather than
%   in whatever arbitrary order they were originally listed).
%
% Inputs:
%   x, y  - coordinate vectors of the polygon's vertices (any order)
%
% Output:
%   order - indices that sort the vertices by angle about the centroid

% Centroid of the point set, used as the reference for angular sorting.
cx = mean(x);
cy = mean(y);

% Angle of each point relative to the centroid.
a = atan2(y - cy, x - cx);

% Sorting by angle arranges the points sequentially around the polygon.
[~, order] = sort(a);
end
