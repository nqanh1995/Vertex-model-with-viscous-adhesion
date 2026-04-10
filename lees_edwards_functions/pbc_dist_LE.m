function [dx,dy] = pbc_dist_LE(dx_old,dy_old,box_L,strain_value)
% lees- edwards pbc
% http://homepage.univie.ac.at/franz.vesely/simsp/dx/node40.html
cy=round(dy_old/box_L);
cdelx=dx_old-cy*strain_value*box_L;
dx=cdelx-round(cdelx/box_L)*box_L;
dy=dy_old-cy*box_L;



end