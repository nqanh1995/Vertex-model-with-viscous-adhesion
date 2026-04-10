function [V,C,topo,strain,stress_tensor] = unshear_initial_state_qvm...
    (V,C,KP,P0,KA,A0)
% input must be a strain zero state!
topo =  get_adj_mat(V,C);
box_size = sqrt(topo.Ncell);
max_relaxation_steps = 5000;
force_tol = 1e-8;
strain = 0;
[V,mf] = vm_minimize_LE(V,C,KP,P0,KA,A0,max_relaxation_steps,force_tol,strain);
stress_tensor = get_qvm_stress(V,C,topo,KP,P0,KA,A0,strain);
disp(['sigma_xy=',num2str(stress_tensor(1,2))]);
s0 = stress_tensor(1,2);
delta_strain = 1e-3;
t1_threshold = 1e-3;
max_n_try = 20;
n_try = 0;

while n_try < max_n_try

    strain = strain + delta_strain;
    cy = (V(:,2)-box_size/2)/box_size;
    V(:,1) = V(:,1) + cy*delta_strain*box_size;

    V = vm_minimize_LE(V,C,KP,P0,KA,A0,max_relaxation_steps,force_tol,strain);

    ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
    ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
    [ex,ey] = pbc_dist_LE(ex,ey,box_size,strain);
    edge_L = sqrt(ex.^2 + ey.^2);

    n_t1 = 0;
    while any(edge_L < t1_threshold) && n_t1 < max_n_try
        E0_list = find(edge_L< t1_threshold);
        disp(['number of T1s to be done: ',num2str(length(E0_list))])
        for j = 1:length(E0_list)
            E0 = E0_list(j);
            [V,C,topo] = T1_edgeswap_LE(V,C,topo,E0,strain);
            V = position_mod_LE(V,box_size,strain);


        end
        [V,~] = vm_minimize_LE(V,C,KP,P0,KA,A0,round(max_relaxation_steps/4),force_tol,strain);

        ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
        ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
        [ex,ey] = pbc_dist_LE(ex,ey,box_size,strain);
        edge_L = sqrt(ex.^2 + ey.^2);
        n_t1 = n_t1 + 1;
    end
    [V,mf] = vm_minimize_LE(V,C,KP,P0,KA,A0,max_relaxation_steps,force_tol,strain);





    stress_tensor = get_qvm_stress(V,C,topo,KP,P0,KA,A0,strain);

    s1 = stress_tensor(1,2);
    
    slope_val  = (s1-s0)/delta_strain;
    if isnan(slope_val)
        break
    end
    delta_strain = (-s0/slope_val - delta_strain);
    s0 = s1;
    n_try = n_try + 1;
end

disp(['the final strain = ',num2str(strain)]);
disp(['the max residual force = ', num2str(mf)]);
disp(['the ratio of shear to pressure is = ', num2str(stress_tensor(1,2)/0.6*(stress_tensor(1,1)+stress_tensor(2,2)))]);



end