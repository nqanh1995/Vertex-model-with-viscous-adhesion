function [V, max_force] = vm_minimize_LE(V, C, KP, P0, KA, A0, max_relaxation_steps, force_tol, strain)
% VM_MINIMIZE_LE  Relax vertex positions to (local) mechanical equilibrium
% using the FIRE algorithm, under Lees-Edwards (sheared) periodic
% boundary conditions.
%
%   [V, max_force] = VM_MINIMIZE_LE(V, C, KP, P0, KA, A0,
%   max_relaxation_steps, force_tol, strain) iteratively moves vertices
%   along the net vertex force (from cell area/perimeter elasticity)
%   using FIRE's adaptive time step and mixed velocity/force-direction
%   update, until either the maximum vertex force drops below force_tol
%   or max_relaxation_steps is reached. Positions are kept inside the
%   box using the sheared (Lees-Edwards) wrapping convention.
%
% Inputs:
%   V, C                 - vertex positions and cell-vertex lists
%   KP, P0               - perimeter elasticity coefficient and target
%                          perimeter, per cell
%   KA, A0               - area elasticity coefficient and target area,
%                          per cell
%   max_relaxation_steps - maximum number of FIRE iterations
%   force_tol            - convergence threshold on the max vertex force
%   strain               - current Lees-Edwards shear strain
%
% Outputs:
%   V         - relaxed vertex positions
%   max_force - final maximum vertex force magnitude reached

Ncell = size(C, 1);
box_size = sqrt(Ncell);
Nvert = size(V, 1);

% --- FIRE algorithm parameters (tuned for the sheared/relaxation case;
% slightly different from the non-sheared vm_minimize.m).
Nmin = 5;              % min consecutive "good" steps before accelerating
finc = 1.2;             % time-step growth factor
fdec = 0.3;             % time-step shrink factor on a bad step
f_alpha = 0.99;         % decay factor for the mixing parameter alpha
alpha_start = 0.1;      % initial velocity/force mixing parameter
dt_max = 0.3;           % maximum allowed time step
dt_init = 0.05;         % initial time step
fire_mass = 3;          % effective inertia used in the velocity update

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
    [Vertex_Force, Cell_Perim, Cell_Area] = get_VM_force_LE(V, C, KP, P0, KA, A0, strain);

    Vertex_Force_mag = sqrt(sum(Vertex_Force.^2, 2));
    if any(isnan(Vertex_Force_mag))
        disp('state is inconsistent, exiting!');
        break;
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
    % vertex positions, wrapping back into the box under sheared PBC.
    vel = vel + (dt/fire_mass) * Vertex_Force;
    V = V + dt*vel;
    V = position_mod_LE(V, box_size, strain);

    max_force = max(Vertex_Force_mag);

    % Periodic progress report.
    if mod(i_relax, round(max_relaxation_steps/4)) == 0
        mean_en = mean(KA.*(Cell_Area-A0).^2 + KP.*(Cell_Perim-P0).^2);
        disp(['step = ', num2str(i_relax)]);
        disp(['mean energy = ', num2str(mean_en)]);
        disp(['max res. force  ', num2str(max_force)]);
    end

    i_relax = i_relax + 1;
end
disp(['final max res. force  ', num2str(max_force)]);
end
