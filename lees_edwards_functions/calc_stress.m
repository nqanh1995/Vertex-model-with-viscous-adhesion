function stress_tensor = calc_stress(V,C,topo,tensions,cell_pressures,strain)

% general calculation of stress tensor given the tension and pressures
[cell_perimeters, cell_areas] = get_cell_shape_LE(V,C,strain);

box_size = sqrt(topo.Ncell);

ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
[ex,ey] = pbc_dist_LE(ex,ey,box_size,strain);
edge_L = sqrt(ex.^2 + ey.^2);


% compute tissue-level stress
stress_tensor = zeros(2,2);

% pressure contributions
stress_tensor(1,1) = -sum(cell_pressures.*cell_areas);
stress_tensor(2,2) = -sum(cell_pressures.*cell_areas);

% tension contributions
for i = 1:topo.Nedge
    stress_tensor(1,1) = stress_tensor(1,1) + tensions(i)*ex(i)*ex(i)/edge_L(i);
    stress_tensor(2,2) = stress_tensor(2,2) + tensions(i)*ey(i)*ey(i)/edge_L(i);
    stress_tensor(1,2) = stress_tensor(1,2) + tensions(i)*ex(i)*ey(i)/edge_L(i);
    stress_tensor(2,1) = stress_tensor(2,1) + tensions(i)*ey(i)*ex(i)/edge_L(i);

end
% scaling to make intensive
stress_tensor = stress_tensor / sum(cell_areas);

end

