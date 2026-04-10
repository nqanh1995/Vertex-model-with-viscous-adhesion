function [Vertex_Force, Cell_Perim, Cell_Area] = get_VM_force_LE(V,C,KP,P0,KA,A0,strain)

V_count = size(V,1);
N = length(C);
box_size = sqrt(N);

[Cell_Perim, Cell_Area] = get_cell_shape_LE(V,C,strain);

% partial derivative of cell energies with respect to vertex positions;
dE_H_mat = zeros(N,2*V_count);
for i = 1:N % loop over cells
    neibs = C{i};
    z = length(neibs);
    Vx = V(neibs,1);
    Vy = V(neibs,2);
    % make sure all vertices are within the same PBC image
    [dx,dy] = pbc_dist_LE(Vx-Vx(end),Vy-Vy(end),box_size,strain);
    Vx = Vx(end) + dx;
    Vy = Vy(end) + dy;

            

    
    % make the neigbor list cyclic
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];
    for n = 2:z+1
        Ehx = KA(i)*(Cell_Area(i)-A0(i))*(Vy(n+1)-Vy(n-1));
        Ehy = KA(i)*(Cell_Area(i)-A0(i))*(Vx(n-1)-Vx(n+1));
        minus_x = 2*KP(i)*(Cell_Perim(i)-P0(i))*...
             (Vx(n)-Vx(n-1))/sqrt((Vx(n-1)-Vx(n))^2+(Vy(n-1)-Vy(n))^2);
        plus_x = 2*KP(i)*(Cell_Perim(i)-P0(i))*...
             (Vx(n)-Vx(n+1))/sqrt((Vx(n+1)-Vx(n))^2+(Vy(n+1)-Vy(n))^2);
        
        minus_y = 2*KP(i)*(Cell_Perim(i)-P0(i))*...
             (Vy(n)-Vy(n-1))/sqrt((Vx(n-1)-Vx(n))^2+(Vy(n-1)-Vy(n))^2);
        plus_y = 2*KP(i)*(Cell_Perim(i)-P0(i))*...
             (Vy(n)-Vy(n+1))/sqrt((Vx(n+1)-Vx(n))^2+(Vy(n+1)-Vy(n))^2);

        if ~isnan(minus_x) && ~isnan(minus_y)
            Ehx = Ehx + minus_x;
            Ehy = Ehy + minus_y;
        end
        if ~isnan(plus_x) && ~isnan(plus_y)
            Ehx = Ehx + plus_x;
            Ehy = Ehy + plus_y;
        end
             
       dE_H_mat(i,neibs(n)) = Ehx;
       dE_H_mat(i,neibs(n)+V_count) = Ehy;
    end
end
Vertex_Force = -sum(dE_H_mat,1);
Vertex_Force = [Vertex_Force(1:V_count);Vertex_Force(V_count+1:end)]';

end