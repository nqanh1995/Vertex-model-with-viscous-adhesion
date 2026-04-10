function [new_x, new_y] = drifting_correction(x,y)
%% Correcting for Global Drift
global_dx = diff(x,1,1);
global_dy = diff(y,1,1);
mean_global_dx = mean(global_dx,2);
mean_global_dy = mean(global_dy,2);
non_drift_dx = bsxfun(@minus,global_dx,mean_global_dx);
non_drift_dy = bsxfun(@minus,global_dy,mean_global_dy);

new_x = bsxfun(@plus,cumsum(non_drift_dx),x(1,:));
new_y = bsxfun(@plus,cumsum(non_drift_dy),y(1,:));
end