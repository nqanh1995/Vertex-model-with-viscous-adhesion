function harmonic_oscillatory_shear(strain_amp,nsteps, interval, seed, Npts, P0_value, KA_value, dt, ...
L_thresh, v0_value, Dr_value, omega, B,stress_only)
    

    % Add required functions to the path
    addpath('matlab_functions/');
    addpath('lees_edwards_functions/');

    

     %%%% Set random seed for reproducibility
    s = RandStream('mcg16807','Seed',seed);
    RandStream.setGlobalStream(s);

    box_size = sqrt(Npts);
    cx = rand(Npts,1);
    cy = rand(Npts,1);
    cx = box_size*cx;
    cy = box_size*cy;
    xy=[cx cy];
    [V, C] = make_voronoi_pbc_fast(xy(:,1)/box_size, xy(:,2)/box_size);

    topo = get_adj_mat(V, C);
    KP = ones(topo.Ncell, 1);
    KA = KA_value * ones(topo.Ncell, 1);
    A0 = ones(topo.Ncell, 1);
    P0 = P0_value * ones(topo.Ncell, 1);
    v0_list = v0_value * ones(topo.Nvert, 1);
    Dr_list = Dr_value * ones(topo.Ncell, 1);
    theta = 2 * pi * rand(topo.Ncell, 1);

    
   
    % Data structure for saving results
    data = struct();
    data.parameters = struct();
    data.parameters.nsteps = nsteps;
    data.parameters.interval = interval;
    data.parameters.seed = seed;
    data.parameters.Npts = Npts;
    data.parameters.P0_value = P0_value;
    data.parameters.KA_value = KA_value;
    data.parameters.dt = dt;
    data.parameters.L_thresh = L_thresh;
    data.parameters.v0_value = v0_value;
    data.parameters.Dr_value = Dr_value;
    data.new_vm_info = struct('newV', [], 'newC', []);

 
    total_strain = 0;
    if stress_only == 1
        output_file_name = sprintf('stress_only_B%g_sA%g_Npts%i_seed%d_p%g_KA%g_dt%g_v%g_Dr%g_gamdot%g.mat', ...
                               B,strain_amp, Npts, seed, P0_value, KA_value, dt, v0_value, Dr_value,omega);
    else
        output_file_name = sprintf('B%g_sA%g_Npts%i_seed%d_p%g_KA%g_dt%g_v%g_Dr%g_gamdot%g.mat', ...
                               B,strain_amp ,Npts, seed, P0_value, KA_value, dt, v0_value, Dr_value,omega);
    end

    %%%% obtain initial config, either by using known state or relaxed from
    %%%% random state
        output_step=1;
        strain = 0;
        t=0;
        max_relaxation_steps=10000;
        force_tol = 10^-9;
        [V,~] = vm_minimize_LE(V,C,KP,P0,KA,A0,max_relaxation_steps,force_tol,strain);
        ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
        ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
        [ex,ey] = pbc_dist_LE(ex,ey,box_size,strain);
        edge_L = sqrt(ex.^2 + ey.^2);
    
        n_t1 = 0;
        while any(edge_L < L_thresh) && n_t1 < 500
            E0_list = find(edge_L< L_thresh);
    %         disp(['number of T1s to be done: ',num2str(length(E0_list))])
            for j = 1:length(E0_list)
                E0 = E0_list(j);
                [V,C,topo] = T1_edgeswap_LE(V,C,topo,E0,strain);
                V = position_mod_LE(V,box_size,strain);
                n_t1 = n_t1 + 1;
    
            end
            [V,~] = vm_minimize_LE(V,C,KP,P0,KA,A0,round(max_relaxation_steps),force_tol,strain);
    
            ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
            ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
            [ex,ey] = pbc_dist_LE(ex,ey,box_size,strain);
            edge_L = sqrt(ex.^2 + ey.^2);
            n_t1 = n_t1 + 1;
        end
        [V,~] = vm_minimize_LE(V,C,KP,P0,KA,A0,max_relaxation_steps,10^-9,0);
        cyclic_switch = 1;
   

tic;



    for n = 1:nsteps
        if total_strain>=strain_amp*4*30
            break
        end
        old_strain = strain;
        % apply shear
        cy = (V(:,2)-box_size/2)/box_size;
        strain = strain_amp*sin(omega*(t+dt));
        t= t+dt;
        delta_gamma = strain-old_strain;
        shearing = cy*delta_gamma*box_size;
        V(:,1) = V(:,1) + shearing;
        V = position_mod_LE(V,box_size,mod(strain,1));

        % get forces
        [VF, ~, Cell_Area] = get_VM_force_LE(V, C, KP, P0, KA, A0,mod(strain,1));




        %%%%%%%%%%%%%%%%%%%%%% Check if sum of Cell_Area exceeds the box size squared with precision
        precision = 4; % Number of decimal places to consider
        if round(sum(Cell_Area), precision) > round(box_size^2, precision)
            disp(['Stopping simulation: Sum of Cell_Area (', num2str(sum(Cell_Area)), '), exceeds the box size squared (', num2str(box_size^2), ').']);
            
            data.last_step_newV_before = V;
            data.last_step_T1_newC_before = C;
            data.last_step_strain = strain;
            data.last_step_T1_timer_before=T1_timer;
            draw_edges_LE(V,topo,strain,1);
            draw_state_LE(V,C,strain,2);
            total_strain = 30;



                    % Polarization Dynamics (VECTORIZED)
                theta = theta + randn(topo.Ncell, 1) .* sqrt(2 * dt .* Dr_list);
                [theta_v, ~] = get_theta_vertex(topo, V, theta,mod(strain,1));
        
                % Cell Dynamics
                VF = VF + v0_list .* theta_v;
                V(:,1) = V(:,1) + VF(:,1)*dt;
                V(:,2) = V(:,2) + VF(:,2)*dt;
                V = position_mod_LE(V,box_size,mod(strain,1));
     
             
            data.last_step_newV_after = V;
            data.last_step_T1_newC_after = C;
            data.last_step_T1_timer_after=T1_timer;
            draw_edges_LE(V,topo,strain,3);
            draw_state_LE(V,C,strain,4);
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Polarization Dynamics (VECTORIZED)
        theta = theta + randn(topo.Ncell, 1) .* sqrt(2 * dt .* Dr_list);
        [theta_v, ~] = get_theta_vertex(topo, V, theta,mod(strain,1));

        % Cell Dynamics
        VF = VF + v0_list .* theta_v;
        [vvx,vvy] = get_vertex_velo_viscous(V,topo,VF,B,strain);
        Fvx = vvx - VF(:,1);
        Fvy = vvy - VF(:,2);
        stress = get_viscous_qvm_stress(V,C,topo,KP,P0,KA,A0,strain,Fvx,Fvy);

        
      


        %%%%%%%%%% save data
        
        
    
        % Save vertex positions for each cell, if mod(n, interval) == 0
        if stress_only == 1
            
            output_step = output_step + 1;
            data.strain(output_step)=strain;
            data.sxy(output_step)=stress(1,2);
            data.sxx(output_step)=stress(1,1);
            data.syy(output_step)=stress(2,2);
            data.time(output_step)=t;
            data.last_newV = V;
            data.last_newC = C;
            data.last_theta =theta;
            data.last_cyclic_switch = cyclic_switch ;
            data.total_strain = total_strain;
            data.last_t = t;
        else
            if mod(n, interval) == 0
                output_step = output_step + 1;
                data.new_vm_info(output_step).newV = V;
                data.new_vm_info(output_step).newC = C;
                data.new_vm_info(output_step).strain = strain;
                data.new_vm_info(output_step).sxy = stress(1,2);
                data.new_vm_info(output_step).sxx = stress(1,1);
                data.new_vm_info(output_step).syy = stress(2,2);
                data.last_T1_timer =T1_timer;
                data.last_theta =theta;
                data.last_cyclic_switch = cyclic_switch ;
                data.total_strain=total_strain;
                data.last_strain = strain;
                
            end
        end

    V(:,1) = V(:,1) + vvx*dt;
    V(:,2) = V(:,2) + vvy*dt;
    total_strain = total_strain + abs(delta_gamma);
    V = position_mod_LE(V,box_size,mod(strain,1));

    
    end
    elapsed_time = toc;

    % Create output file name
 
    
    data.Elapsed_time = elapsed_time;
    save(output_file_name, 'data', '-v7.3');
    format long
    fprintf('strain %.9f', total_strain)
end
