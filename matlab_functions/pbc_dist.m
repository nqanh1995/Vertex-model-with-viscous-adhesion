function new_dist = pbc_dist(old_dist,box_size)
pbc_switch = abs(old_dist)>0.5*box_size;
new_dist = (1-pbc_switch).*old_dist + pbc_switch.*(old_dist-sign(old_dist).*box_size);
return