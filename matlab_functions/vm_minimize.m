function V =  vm_minimize(V,C,KP,P0,KA,A0,max_relaxation_steps,force_tol)

Ncell = size(C,1);
box_size = sqrt(Ncell);
Nvert = size(V,1);





Nmin = 5;
finc = 1.1;
fdec = 0.5;
f_alpha = 0.99;
alpha_start = 0.1;
dt_max = 100;
dt_init = 0.01;
fire_mass = 4; % inertia


tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % BEGIN FIRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = dt_init;
vel = zeros(Nvert,2);

alpha = alpha_start;

max_force = 1e8;
i_relax = 0;
Nneg = 0;

while max_force > force_tol && i_relax < max_relaxation_steps
     [Vertex_Force, Cell_Perim, Cell_Area] = get_VM_force(V,C,KP,P0,KA,A0);

    Vertex_Force_mag = sqrt(sum(Vertex_Force.^2,2));
    if any(isnan(Vertex_Force_mag))
        disp('state is inconsistent, exiting!');
        return;
    end
    Vertex_Force_normalized = bsxfun(@rdivide, Vertex_Force, Vertex_Force_mag);



    Pfire = sum(dot(vel,Vertex_Force));





    if Pfire > 0 && Nneg > Nmin
        dt = min(dt*finc,dt_max);
        alpha = alpha*f_alpha;
        vmag = vecnorm(vel,2,2);
        vel = (1-alpha)*vel + alpha* Vertex_Force_normalized .* vmag;
    end
    if Pfire <= 0
        dt = dt*fdec;
        alpha = alpha_start;
        vel = zeros(Nvert,2);
        Nneg = 0;
    end
    Nneg = Nneg+1;
    vel = vel + (dt/fire_mass)*Vertex_Force;
    V = mod(V + dt*vel,box_size);

    max_force = max(Vertex_Force_mag);

    if mod(i_relax,round(max_relaxation_steps/10)) == 0
        mean_en = mean(KA.*(Cell_Area-A0).^2 + KP.*(Cell_Perim-P0).^2);
        disp(['step = ', num2str(i_relax)]);
        disp(['mean energy = ', num2str(mean_en)]);
        disp(['mean force  ', num2str(max_force)]);
        % disp(['min sa / mean sa= ', num2str(min(cell_sa)/mean(cell_sa))]);
        % disp(['current dt=',num2str(dt)]);
        % disp(['std volume / mean volume= ', num2str(std(cell_vol)/V0(1))]);


    end


    i_relax = i_relax + 1;
end

end
