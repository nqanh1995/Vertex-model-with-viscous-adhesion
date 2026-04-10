function [theta_v, edge_L] = get_theta_vertex(topo, V, theta,strain)

box_size = sqrt(topo.Ncell);
Ncell = topo.Ncell;
Nvert = topo.Nvert;

edge_dx = V(topo.E_list(:,1),1) - V(topo.E_list(:,2),1);
edge_dy = V(topo.E_list(:,1),2) - V(topo.E_list(:,2),2);
[edge_dx,edge_dy]=pbc_dist_LE(edge_dx,edge_dy,box_size,strain);
edge_L = sqrt(edge_dx.^2 + edge_dy.^2);
L = topo.VEadj*edge_L;

theta_c = [cos(theta), sin(theta)];

theta_v = zeros(Nvert,2);
for i = 1:Nvert
        neibs_C = find(topo.CVadj(:,i));
        z = length(neibs_C);
        n_c = zeros(z,2);
        n_c(:,1) = theta_c(neibs_C,1);
        n_c(:,2) = theta_c(neibs_C,2);
        a_c = zeros(z,1);
        E_Vi = find(topo.VEadj(i,:));  % Find common edges between cells
     for j = 1:z
            cid = neibs_C(j);
            E_Cj = find(topo.CEadj(cid,:));
            l12 = intersect(E_Cj,E_Vi);
            a_c(j)=sum(edge_L(l12))./(2*z*L(i));
     end

         theta_v(i,:)=sum(a_c.*n_c, 1);

        % Normalize theta_v to make it a unit vector
        magnitude = norm(theta_v(i,:));
        if magnitude > 0
            theta_v(i,:) = theta_v(i,:) / magnitude;
        else
            warning('theta_v at vertex %d has zero magnitude and cannot be normalized.', i);
        end
end

end



