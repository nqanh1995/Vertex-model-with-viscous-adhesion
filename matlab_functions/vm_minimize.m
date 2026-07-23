function V = vm_minimize(V, C, KP, P0, KA, A0, max_relaxation_steps, force_tol)
% VM_MINIMIZE  Relax vertex positions to (local) mechanical equilibrium
% using the FIRE (Fast Inertial Relaxation Engine) algorithm.
%
%   V = VM_MINIMIZE(V, C, KP, P0, KA, A0, max_relaxation_steps, force_tol)
%   iteratively moves vertices along the net vertex force (from cell
%   area/perimeter elasticity) using FIRE's adaptive time step and mixed
%   velocity/force-direction update, until either the maximum vertex
%   force drops below force_tol or max_relaxation_steps is reached.
%
% Inputs:
%   V, C                 - vertex positions and cell-vertex lists
%   KP, P0               - perimeter elasticity coefficient and target
%                          perimeter, per cell
%   KA, A0               - area elasticity coefficient and target area,
%                          per cell
%   max_relaxation_steps - maximum number of FIRE iterations
%   force_tol            - convergence threshold on the max vertex force
%
% Output:
%   V - relaxed vertex positions
%
% Note: this uses GET_VM_FORCE_LE (the Lees-Edwards / sheared-PBC force
% calculation) evaluated at strain = 0, which is exactly equivalent to a
% flat (non-sheared) periodic box; this avoids depending on a separate
% non-sheared GET_VM_FORCE function. Requires 'lees_edwards_functions/'
% to be on the path.

Ncell = size(C, 1);
box_size = sqrt(Ncell);
Nvert = size(V, 1);

% --- FIRE algorithm parameters (standard values).
Nmin = 5;              % min consecutive "good" steps before accelerating
finc = 1.1;             % time-step growth factor
fdec = 0.5;             % time-step shrink factor on a bad step
f_alpha = 0.99;         % decay factor for the mixing parameter alpha
alpha_start = 0.1;      % initial velocity/force mixing parameter
dt_max = 1;             % maximum allowed time step
dt_init = 0.01;         % initial time step
fire_mass = 4;          % effective inertia used in the velocity update

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % BEGIN FIRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = dt_init;
vel = zeros(Nvert, 2);
alpha = alpha_start;

max_force = 1e8;
i_relax = 0;
Nneg = 0;   % consecutive number of "good" (power > 0) steps

while max_force > force_tol && i_relax < max_relaxation_steps
    % Net elastic force on every vertex, plus current cell shape info.
    % (get_VM_force_LE at strain = 0 is exactly the flat/non-sheared
    % force calculation.)
    [Vertex_Force, Cell_Perim, Cell_Area] = get_VM_force_LE(V, C, KP, P0, KA, A0, 0);

    Vertex_Force_mag = sqrt(sum(Vertex_Force.^2, 2));
    if any(isnan(Vertex_Force_mag))
        disp('state is inconsistent, exiting!');
        return;
    end
    Vertex_Force_normalized = bsxfun(@rdivide, Vertex_Force, Vertex_Force_mag);

    % FIRE "power": how aligned the current velocity is with the force.
    Pfire = sum(dot(vel, Vertex_Force));

    if Pfire > 0 && Nneg > Nmin
        % Velocity is well-aligned with the force: speed up (larger dt)
        % and bias the velocity further toward the force direction.
        dt = min(dt*finc, dt_max);
        alpha = alpha*f_alpha;
        vmag = vecnorm(vel, 2, 2);
        vel = (1-alpha)*vel + alpha*Vertex_Force_normalized .* vmag;
    end
    if Pfire <= 0
        % Overshot / moving against the force: reset velocity and time
        % step to restart the descent more conservatively.
        dt = dt*fdec;
        alpha = alpha_start;
        vel = zeros(Nvert, 2);
        Nneg = 0;
    end
    Nneg = Nneg + 1;

    % Semi-implicit Euler update: accelerate by the force, then advance
    % vertex positions, wrapping back into the periodic box.
    vel = vel + (dt/fire_mass) * Vertex_Force;
    V = mod(V + dt*vel, box_size);

    max_force = max(Vertex_Force_mag);

    % Periodic progress report.
    if mod(i_relax, round(max_relaxation_steps/10)) == 0
        mean_en = mean(KA.*(Cell_Area-A0).^2 + KP.*(Cell_Perim-P0).^2);
        disp(['step = ', num2str(i_relax)]);
        disp(['mean energy = ', num2str(mean_en)]);
        disp(['mean force  ', num2str(max_force)]);
    end

    i_relax = i_relax + 1;
end
end
