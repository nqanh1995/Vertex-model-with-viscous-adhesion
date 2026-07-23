function avm_viscous(nsteps, interval, seed, Npts, P0_value, KA_value, dt, tau_T1, L_thresh, v0_value, Dr_value, B)
% AVM_VISCOUS  Run an Active Vertex Model (AVM) simulation of a confluent
% epithelial tissue with cell-cell viscous adhesion, on a flat (non-
% sheared) periodic box, and save the results to a .mat file.
%
%   AVM_VISCOUS(nsteps, interval, seed, Npts, P0_value, KA_value, dt,
%   tau_T1, L_thresh, v0_value, Dr_value, B) builds a random Voronoi
%   tissue of Npts cells, relaxes it to mechanical equilibrium, then
%   advances the simulation for nsteps time steps under active
%   self-propulsion and viscous vertex-vertex coupling, performing T1
%   (neighbor-exchange) transitions as edges shrink below threshold.
%   Snapshots are saved every "interval" steps.
%
% Inputs:
%   nsteps     - number of simulation steps to run
%   interval   - how often (in steps) to save a snapshot of the tissue
%   seed       - random seed, for reproducibility
%   Npts       - number of cells in the system
%   P0_value   - target shape index p0, applied uniformly to all cells
%   KA_value   - area elasticity coefficient, applied uniformly
%   dt         - integration time step
%   tau_T1     - cooldown time an edge must wait between T1 transitions
%   L_thresh   - edge length threshold that triggers a T1 transition
%   v0_value   - active self-propulsion speed, applied uniformly
%   Dr_value   - rotational diffusion constant of the polarity angle,
%                applied uniformly
%   B          - cell-cell viscous coefficient (xi)
%
% Output: none (results are saved to a .mat file named after the
% parameters, in the current working directory).

% Add required functions to the path.
addpath('matlab_functions/');
addpath('lees_edwards_functions/');

% Set random seed for reproducibility.
s = RandStream('mcg16807', 'Seed', seed);
RandStream.setGlobalStream(s);

%% Build the initial tissue geometry: random points -> periodic Voronoi.
box_size = sqrt(Npts);
cx = rand(Npts, 1);
cy = rand(Npts, 1);
cx = box_size * cx;
cy = box_size * cy;
xy = [cx cy];
[V, C] = make_voronoi_pbc_fast(xy(:,1)/box_size, xy(:,2)/box_size);

%% Initial parameters (uniform across all cells/vertices).
topo = get_adj_mat(V, C);
KP = ones(topo.Ncell, 1);
KA = KA_value * ones(topo.Ncell, 1);
A0 = ones(topo.Ncell, 1);
P0 = P0_value * ones(topo.Ncell, 1);
v0_list = v0_value * ones(topo.Nvert, 1);
Dr_list = Dr_value * ones(topo.Ncell, 1);
theta = 2 * pi * rand(topo.Ncell, 1);   % initial cell polarity angles

% T1 timers (time since each edge's last T1 transition) and bookkeeping
% for observed waiting times / event counts.
T1_timer = rand(topo.Nedge, 1);
tau_wait = zeros(topo.Nedge, 1);
nT1_edge = zeros(topo.Nedge, 1);
nT1_events = zeros(nsteps, 1);
all_tau_wait = [];  % Collect all observed waiting times

%% Data structure for saving results.
data = struct();
data.parameters = struct();
data.parameters.nsteps = nsteps;
data.parameters.B = B;
data.parameters.interval = interval;
data.parameters.seed = seed;
data.parameters.Npts = Npts;
data.parameters.P0_value = P0_value;
data.parameters.KA_value = KA_value;
data.parameters.dt = dt;
data.parameters.tau_T1 = tau_T1;
data.parameters.L_thresh = L_thresh;
data.parameters.v0_value = v0_value;
data.parameters.Dr_value = Dr_value;
data.new_vm_info = struct('V', [], 'C', []);

strain = 0;   % this script never applies shear; strain stays at 0 throughout

output_file_name = sprintf('B%g_Npts%i_seed%d_p%g_KA%g_dt%g_v%g_Dr%g_T1_%g.mat', ...
                           B, Npts, seed, P0_value, KA_value, dt, v0_value, Dr_value, tau_T1);

tic;   % start the simulation-wall-clock timer (read out via toc at the end)

%% Initial relaxation: bring the random Voronoi tessellation to mechanical
% equilibrium, resolving any too-short edges with T1 transitions along
% the way, before starting the dynamics.
output_step = 1;
max_relaxation_steps = 5000;
force_tol = 10^-9;
[V, ~] = vm_minimize_LE(V, C, KP, P0, KA, A0, max_relaxation_steps, force_tol, strain);
ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
[ex, ey] = pbc_dist_LE(ex, ey, box_size, strain);
edge_L = sqrt(ex.^2 + ey.^2);

n_t1 = 0;
while any(edge_L < L_thresh) && n_t1 < 500
    E0_list = find(edge_L < L_thresh);
    for j = 1:length(E0_list)
        E0 = E0_list(j);
        [V, C, topo] = T1_edgeswap_LE(V, C, topo, E0, strain);
        V = position_mod_LE(V, box_size, strain);
    end
    [V, ~] = vm_minimize_LE(V, C, KP, P0, KA, A0, round(max_relaxation_steps), force_tol, strain);

    ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
    ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
    [ex, ey] = pbc_dist_LE(ex, ey, box_size, strain);
    edge_L = sqrt(ex.^2 + ey.^2);
    n_t1 = n_t1 + 1;
end
[V, ~] = vm_minimize_LE(V, C, KP, P0, KA, A0, max_relaxation_steps, 10^-9, 0);

% Save the relaxed initial state as the first snapshot.
data.new_vm_info(output_step).V = V;
data.new_vm_info(output_step).C = C;

%% Main simulation loop.
for n = 1:nsteps

    V = position_mod_LE(V, box_size, mod(strain, 1));
    ex = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
    ey = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
    [ex, ey] = pbc_dist_LE(ex, ey, box_size, strain);
    edge_L = sqrt(ex.^2 + ey.^2);

    %%%%%% T1 transition: flip any edge that has both shrunk below
    % L_thresh and waited at least tau_T1 since its last transition.
    E0_list = find(edge_L < L_thresh);
    if ~isempty(E0_list)
        for jj = 1:length(E0_list)
            E0 = E0_list(jj);
            if T1_timer(E0) > tau_T1
                [V, C, topo] = T1_edgeswap_LE(V, C, topo, E0, mod(strain, 1));

                % Update waiting times and reset timer for involved cells.
                tau_wait(E0) = tau_wait(E0) + T1_timer(E0);
                nT1_edge(E0) = nT1_edge(E0) + 1;
                T1_timer(E0) = 0;

                % Record T1 event.
                nT1_events(n) = nT1_events(n) + 1;

                disp(['Performed T1 on edge ', num2str(E0), ', with length of ', num2str(edge_L(E0))]);
            end
        end
    end

    % Get mechanical (area + perimeter elasticity) forces on each vertex.
    [VF, ~, Cell_Area] = get_VM_force_LE(V, C, KP, P0, KA, A0, mod(strain, 1));

    %%%%%%%%%%%%%%%%%%%%%% Sanity check: total cell area should never
    % exceed the box area. If it does, something has gone numerically
    % wrong (e.g. an inverted/degenerate cell), so the simulation state
    % just before and after one more explicit (non-viscous) update step
    % is saved for debugging, along with diagnostic plots, and the run
    % is stopped early.
    precision = 4; % Number of decimal places to consider
    if round(sum(Cell_Area), precision) > round(box_size^2, precision)
        disp(['Stopping simulation: Sum of Cell_Area (', num2str(sum(Cell_Area)), '), exceeds the box size squared (', num2str(box_size^2), ').']);

        data.last_step_newV_before = V;
        data.last_step_T1_newC_before = C;
        data.last_step_strain = strain;
        data.last_step_T1_timer_before = T1_timer;
        draw_edges_LE(V, topo, strain, 1);    %%% plot tissue state
        draw_state_LE(V, C, strain, 2);       %%% plot tissue state

        % Polarization Dynamics (VECTORIZED)
        theta = theta + randn(topo.Ncell, 1) .* sqrt(2 * dt .* Dr_list);
        [theta_v, ~] = get_theta_vertex(topo, V, theta, mod(strain, 1));

        % Cell Dynamics (simple explicit Euler step, no viscous coupling,
        % just to record the state one step past the failure point).
        VF = VF + v0_list .* theta_v;
        V(:,1) = V(:,1) + VF(:,1)*dt;
        V(:,2) = V(:,2) + VF(:,2)*dt;
        V = position_mod_LE(V, box_size, mod(strain, 1));

        % Update T1 timers.
        T1_timer = T1_timer + dt;

        data.last_step_newV_after = V;
        data.last_step_T1_newC_after = C;
        data.last_step_T1_timer_after = T1_timer;
        draw_edges_LE(V, topo, strain, 3);       %%% plot tissue state
        draw_state_LE(V, C, strain, 4);          %%% plot tissue state
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Polarization Dynamics (VECTORIZED): each cell's polarity angle
    % undergoes rotational diffusion, then is projected onto vertices.
    theta = theta + randn(topo.Ncell, 1) .* sqrt(2 * dt .* Dr_list);
    [theta_v, ~] = get_theta_vertex(topo, V, theta, mod(strain, 1));

    % Cell Dynamics: add the active self-propulsion force to the
    % mechanical force, then solve for vertex velocities including the
    % cell-cell viscous coupling.
    VF = VF + v0_list .* theta_v;
    [vvx, vvy] = get_vertex_velo_viscous(V, topo, VF, B, strain);

    %%%%%%%%%% Save data: record a snapshot every "interval" steps.
    if mod(n, interval) == 0
        data.new_vm_info(output_step).V = V;
        data.new_vm_info(output_step).C = C;
        data.new_vm_info(output_step).VF = VF;
        output_step = output_step + 1;
    end

    % Advance vertex positions by one explicit Euler step.
    V(:,1) = V(:,1) + vvx*dt;
    V(:,2) = V(:,2) + vvy*dt;
    V = position_mod_LE(V, box_size, mod(strain, 1));

    % Update T1 timers.
    T1_timer = T1_timer + dt;
end

elapsed_time = toc;

data.Elapsed_time = elapsed_time;
save(output_file_name, 'data', '-v7.3');
end
