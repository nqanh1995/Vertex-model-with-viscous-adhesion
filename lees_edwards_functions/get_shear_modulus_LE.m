function [sm,dE2_mat] = get_shear_modulus_LE(Npts,P0,V,newC,gamma,ka)




V_count = size(V,1); 
box_size = sqrt(Npts);
P0 = P0*ones(Npts,1);
K_A = ka*ones(Npts,1); % set modulus for area term to be zero
K_P = ones(Npts,1); % set modulus for area term to be one
A0 = ones(Npts,1); 

gamma = mod(gamma,1);


dE2_mat = get_dE2_mat_vm_LE(V,newC,K_A,K_P,A0,P0,gamma);
% [pe, ar] = get_cell_shape(V,newC);
[pe, ar] = get_cell_shape_LE(V,newC,gamma);
affine_part = zeros(Npts,1);
dE2_gamma_V = zeros(Npts,2*V_count);
zlist = zeros(Npts,1);

for i = 1:Npts % loop over cells
    neibs = newC{i};
    z = length(neibs);
    zlist(i) = z;
    Vx = V(neibs,1);
    Vy = V(neibs,2);
    dx = Vx(2:z)-Vx(1);
    dy = Vy(2:z)-Vy(1);
    [dx,dy]=pbc_dist_LE(dx,dy,box_size,gamma);
    Vx(2:z) = Vx(1) + dx;
    Vy(2:z) = Vy(1) + dy;

%     Vx(2:z) = Vx(1) + pbc_dist(Vx(2:z)-Vx(1),box_size);
%     Vy(2:z) = Vy(1) + pbc_dist(Vy(2:z)-Vy(1),box_size);
    % make sure all vertices are within the same PBC image
%     Vx = cx(i) + pbc_dist(V(neibs,1)-cx(i),box_size);
%     Vy = cy(i) + pbc_dist(V(neibs,2)-cy(i),box_size);
%     [~,ord] = sort(mod(atan2(Vy-mean(Vy),Vx-mean(Vx)),2*pi));
%     neibs = neibs(ord);
%     Vx = Vx(ord);
%     Vy = Vy(ord);
%     
    % make the neigbor list cyclic
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];

    part0 = 0;
    part1 = 0;
    part2 = 0;
    LL = zeros(z+2,1); % edge lengths
    for n = 1:z
       LL(n) = sqrt((Vx(n+1)-Vx(n))^2+(Vy(n+1)-Vy(n))^2);
       part0 = part0 + (-Vx(n+1) + Vx(n))*((-Vy(n+1) + Vy(n)))/LL(n);
       part1 = part1 + -(-Vx(n+1) + Vx(n))^2*(-Vy(n+1) + Vy(n))^2/LL(n)^3;
       part2 = part2 + (-Vy(n+1) + Vy(n))^2/LL(n);

    end
    LL(z+1) = LL(1);
    LL(z+2) = LL(2);
    affine_part(i) = 2*K_P(i)*(part0^2 + (part1+part2)*(pe(i)-P0(i)));
    
    
    % calculate partial^2 E_i / (partial gamma partial V)
    
    for n = 2:z+1
        dE2_gamma_V(i,neibs(n)) = ...
            2*K_P(i)*(part0*((Vx(n) - Vx(n+1))/LL(n) + (Vx(n) - ...
            Vx(n-1))/LL(n-1)) + (-P0(i) + pe(i))*((Vy(n) - Vy(n+1))/LL(n) + ...
            ((Vx(n) - Vx(n+1))^2*(-Vy(n) + Vy(n+1)))/LL(n)^3 + (Vy(n) - Vy(n-1))/LL(n-1) - ...
            ((Vx(n) - Vx(n-1))^2*(Vy(n) - Vy(n-1)))/LL(n-1)^3));
        dE2_gamma_V(i,neibs(n)+V_count) = ...
            2*K_P(i)*(part0*((Vy(n) - Vy(n+1))/LL(n) + (Vy(n) - Vy(n-1))/LL(n-1)) + (-P0(i) + pe(i))*((Vx(n) - Vx(n+1))/LL(n) + (Vx(n) - Vx(n-1))/LL(n-1) + ((-Vx(n) + Vx(n+1))*(Vy(n) - Vy(n+1))^2)/LL(n)^3 - ((Vx(n) - Vx(n-1))*(Vy(n) - Vy(n-1))^2)/LL(n-1)^3));  
        
    end

%     for n = 2:z+1
%         dE2_gamma_V(i,neibs(n)) = ...
%             2*K_P(i)*(part0*((Vx(n) - Vx(n+1))/LL(n) + (Vx(n) - ...
%             Vx(n-1))/LL(n-1)) + (-P0(i) + pe(i))*((Vy(n) - Vy(n+1))/LL(n) + ...
%             ((Vx(n) - Vx(n+1))^2*(-Vy(n) + Vy(n+1)))/LL(n)^3 + (Vy(n) - Vy(n-1))/LL(n-1) - ...
%             ((Vx(n) - Vx(n-1))^2*(Vy(n) - Vy(n-1)))/LL(n-1)^3));
%         dE2_gamma_V(i,neibs(n)+2*Npts) = ...
%             2*K_P(i)*(part0*((Vy(n) - Vy(n+1))/LL(n) + (Vy(n) - Vy(n-1))/LL(n-1)) + (-P0(i) + pe(i))*((Vx(n) - Vx(n+1))/LL(n) + (Vx(n) - Vx(n-1))/LL(n-1) + ((-Vx(n) + Vx(n+1))*(Vy(n) - Vy(n+1))^2)/LL(n)^3 - ((Vx(n) - Vx(n-1))*(Vy(n) - Vy(n-1))^2)/LL(n-1)^3));  
        
%     end
        
end


dE2_gamma_V = sum(dE2_gamma_V,1);
shear_mod_affine = sum(affine_part)/Npts;
shear_mod_nonaffine =  (dE2_gamma_V)*(dE2_mat\(dE2_gamma_V)')/Npts;
sm = shear_mod_affine - shear_mod_nonaffine;

end