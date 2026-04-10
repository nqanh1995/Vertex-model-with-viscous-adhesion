function [stress_tensor,cell_stress] = get_viscous_qvm_stress(V,C,topo,KP,P0,KA,A0,strain,Fvx,Fvy)

[cell_perimeters, cell_areas] = get_cell_shape_LE(V,C,strain);
% tensions and pressures
tensions = (2*KP.*(cell_perimeters  - P0))'*topo.CEadj;
tensions = tensions(:);

cell_pressures = -2*KA.*(cell_areas-A0);
cell_pressures = cell_pressures(:);

% stress tensor of each cell, ordered as sxx, sxy, syx, syy
cell_stress = zeros(topo.Ncell,4);

box_size = sqrt(topo.Ncell);

ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
[ex,ey] = pbc_dist_LE(ex,ey,box_size,strain);
edge_L = sqrt(ex.^2 + ey.^2);


% tissue-level stress
stress_tensor = zeros(2,2);

% pressure contributions
for c=1:topo.Ncell
    cell_stress(c,1)=-cell_pressures(c);
    cell_stress(c,4)=-cell_pressures(c);
end

% tension contributions

for c=1:topo.Ncell
    eid = topo.CEadj(c,:);
    cell_stress(c,1)= cell_stress(c,1) + ...
        sum(tensions(eid).*ex(eid).*ex(eid)./edge_L(eid))/(2*cell_areas(c));
    cell_stress(c,2)= cell_stress(c,2) + ...
        sum(tensions(eid).*ex(eid).*ey(eid)./edge_L(eid))/(2*cell_areas(c));
    cell_stress(c,3)= cell_stress(c,3) + ...
        sum(tensions(eid).*ey(eid).*ex(eid)./edge_L(eid))/(2*cell_areas(c));
    cell_stress(c,4)= cell_stress(c,4) + ...
        sum(tensions(eid).*ey(eid).*ey(eid)./edge_L(eid))/(2*cell_areas(c));
end


% viscous force contribution

for c=1:topo.Ncell
    verts = C{c};
    z = length(verts);
    Vx = V(verts,1);
    Vy = V(verts,2);
    [dvx,dvy]=pbc_dist_LE(Vx(2:z)-Vx(1),Vy(2:z)-Vy(1),box_size,strain);
    Vx(2:z) = Vx(1) + dvx;
    Vy(2:z) = Vy(1) + dvy;
    ccx = mean(Vx, 1); 
    ccy = mean(Vy, 1); 
    [Rx,Ry] = pbc_dist_LE(Vx-ccx,Vy-ccy,box_size,strain);
    cell_stress(c,1) = cell_stress(c,1) -sum(Rx.*Fvx(verts) + Rx.*Fvx(verts))/(2*z*cell_areas(c));
    cell_stress(c,2) = cell_stress(c,2) -sum(Rx.*Fvy(verts) + Ry.*Fvx(verts))/(2*z*cell_areas(c));
    cell_stress(c,3) = cell_stress(c,3) -sum(Ry.*Fvx(verts) + Rx.*Fvy(verts))/(2*z*cell_areas(c));
    cell_stress(c,4) = cell_stress(c,4) -sum(Ry.*Fvy(verts) + Ry.*Fvy(verts))/(2*z*cell_areas(c));
end

% adding cell stress to get tissue stress
stress_tensor(1,1) = sum(cell_stress(:,1).*cell_areas);
stress_tensor(1,2) = sum(cell_stress(:,2).*cell_areas);
stress_tensor(2,1) = sum(cell_stress(:,3).*cell_areas);
stress_tensor(2,2) = sum(cell_stress(:,4).*cell_areas);

% scaling to make intensive
stress_tensor = stress_tensor / sum(cell_areas);

end