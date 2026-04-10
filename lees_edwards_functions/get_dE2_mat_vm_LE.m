function dE2_mat = get_dE2_mat_vm_LE(V,newC,K_A,K_P,A0,P0,strain)

% 2nd partial derivative of  energy with respect to vertex positions;
% \frac{\partial^2 E}{\partial h_{m\alpha} \partial h_{n,\beta}}
% see dE2_matrix.nb for analytical calculations

strain=mod(strain,1);
[pe, ar] = get_cell_shape_LE(V,newC,strain);
Npts = length(newC);
box_size = sqrt(Npts);


V_count = size(V,1);

dE2_mat = zeros(2*V_count,2*V_count);

for i = 1:Npts % loop over cells
    dE2_mat_cell = zeros(2*V_count,2*V_count);
    neibs = newC{i};
    z = length(neibs);
    Vx = V(neibs,1);
    Vy = V(neibs,2);

    % make sure all vertices are within the same PBC image 
    dx = Vx(2:z)-Vx(1);
    dy = Vy(2:z)-Vy(1);
    [dx,dy]=pbc_dist_LE(dx,dy,box_size,strain);
    Vx(2:z) = Vx(1) + dx;
    Vy(2:z) = Vy(1) + dy;

%     Vx(2:z) = Vx(1) + pbc_dist(Vx(2:z)-Vx(1),box_size);
%     Vy(2:z) = Vy(1) + pbc_dist(Vy(2:z)-Vy(1),box_size);
%     [~,ord] = sort(mod(atan2(Vy-mean(Vy),Vx-mean(Vx)),2*pi));
%     neibs = neibs(ord);
%     Vx = Vx(ord);
%     Vy = Vy(ord);
    
    % make the neigbor list cyclic
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];


	for m = 2:z+1
		for n = 2:z+1
			if mod(n,z) == mod(m,z)	
				%	m == n
				% x,x
				dE2_mat_cell(neibs(m),neibs(n)) = ...
					2*K_P(i)*(-P0(i) + pe(i))*(-((Vx(m) - Vx(m+1))^2/((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2)^(3/2)) + 1/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) - (Vx(m) - Vx(m-1))^2/((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)^(3/2) + 1/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + ...
	 				2*K_P(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2))^2 + (K_A(i)*(Vy(m+1) - Vy(m-1))^2)/2;

				% y, y
				dE2_mat_cell(neibs(m)+V_count,neibs(n)+V_count) = ...
					(K_A(i)*(Vx(m+1) - Vx(m-1))^2)/2 + 2*K_P(i)*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2))^2 + ...
	 				2*K_P(i)*(-P0(i) + pe(i))*(1/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) - (Vy(m) - Vy(m+1))^2/((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2)^(3/2) + 1/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2) - (Vy(m) - Vy(m-1))^2/((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)^(3/2));

				% x, y == y,x
				dE2_mat_cell(neibs(m),neibs(n)+V_count) = ...
					2*K_P(i)*(-P0(i) + pe(i))*(-(((Vx(m) - Vx(m+1))*(Vy(m) - Vy(m+1)))/((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2)^(3/2)) - ((Vx(m) - Vx(m-1))*(Vy(m) - Vy(m-1)))/((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)^(3/2)) + ...
	 				2*K_P(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2))*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (K_A(i)*(-Vx(m+1) + Vx(m-1))*(Vy(m+1) - Vy(m-1)))/2;

				dE2_mat_cell(neibs(m)+V_count,neibs(n)) = dE2_mat_cell(neibs(m),neibs(n)+V_count);

			elseif mod(n,z) == mod(m+1,z) 
				% n == m+1
				% x,x
				dE2_mat_cell(neibs(m),neibs(n)) = ...
					(-2*K_P(i)*(-P0(i) + pe(i))*(Vy(m) - Vy(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 2*K_P(i)*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	  				((Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (K_A(i)*(-Vy(m) + Vy(n+1))*(Vy(n) - Vy(m-1)))/2;


				% y, y
				dE2_mat_cell(neibs(m)+V_count,neibs(n)+V_count) = ...
					(K_A(i)*(Vx(m) - Vx(n+1))*(-Vx(n) + Vx(m-1)))/2 - (2*K_P(i)*(-P0(i) + pe(i))*(Vx(m) - Vx(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 2*K_P(i)*((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	  				((Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2));

				% x, y 
				dE2_mat_cell(neibs(m),neibs(n)+V_count) = ...
					(-((2*A0(i) - 2*ar(i))*K_A(i)) + (4*K_P(i)*(-P0(i) + pe(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 4*K_P(i)*((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	   				((Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + K_A(i)*(Vx(m) - Vx(n+1))*(Vy(n) - Vy(m-1)))/2;

				% y, x
				dE2_mat_cell(neibs(m)+V_count,neibs(n)) = ...
					((2*A0(i) - 2*ar(i))*K_A(i) + (4*K_P(i)*(-P0(i) + pe(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + K_A(i)*(Vx(n) - Vx(m-1))*(Vy(m) - Vy(n+1)) + 4*K_P(i)*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	   				((Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)))/2;


			elseif mod(n,z) == mod(m-1,z) 
				% x,x
				dE2_mat_cell(neibs(m),neibs(n)) = ...
					2*K_P(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2))*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2)) - ...
	 				(2*K_P(i)*(-P0(i) + pe(i))*(Vy(m) - Vy(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + (K_A(i)*(Vy(m) - Vy(n-1))*(Vy(m+1) - Vy(n)))/2;

				% y, y
				dE2_mat_cell(neibs(m)+V_count,neibs(n)+V_count) = ...
					(K_A(i)*(Vx(m) - Vx(n-1))*(Vx(m+1) - Vx(n)))/2 - (2*K_P(i)*(-P0(i) + pe(i))*(Vx(m) - Vx(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 2*K_P(i)*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2))* ...
					((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2));

				% x, y 
				dE2_mat_cell(neibs(m),neibs(n)+V_count) = ...
					((2*A0(i) - 2*ar(i))*K_A(i) + (4*K_P(i)*(-P0(i) + pe(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + K_A(i)*(-Vx(m) + Vx(n-1))*(Vy(m+1) - Vy(n)) + 4*K_P(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2))* ...
	   				((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2)))/2;

				% y, x
				dE2_mat_cell(neibs(m)+V_count,neibs(n)) = ...
					(-((2*A0(i) - 2*ar(i))*K_A(i)) + K_A(i)*(-Vx(m+1) + Vx(n))*(Vy(m) - Vy(n-1)) + 4*K_P(i)*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2))* ...
	   				((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)) + (4*K_P(i)*(-P0(i) + pe(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2))/2;
				

			else % m and n are not neighboring vertices 
				% x,x
				dE2_mat_cell(neibs(m),neibs(n)) = ...
					2*K_P(i)*((-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (K_A(i)*(-Vy(n-1) + Vy(n+1))*(Vy(m+1) - Vy(m-1)))/2;


				% y, y
				dE2_mat_cell(neibs(m)+V_count,neibs(n)+V_count) = ...
					(K_A(i)*(Vx(n-1) - Vx(n+1))*(-Vx(m+1) + Vx(m-1)))/2 + 2*K_P(i)*((-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2));


				% x, y 
				dE2_mat_cell(neibs(m),neibs(n)+V_count) = ...
					2*K_P(i)*((-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (K_A(i)*(Vx(n-1) - Vx(n+1))*(Vy(m+1) - Vy(m-1)))/2;


				% y, x
				dE2_mat_cell(neibs(m)+V_count,neibs(n)) = ...
					(K_A(i)*(Vx(m+1) - Vx(m-1))*(Vy(n-1) - Vy(n+1)))/2 + 2*K_P(i)*((-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2));

				
			end
		end
	end


dE2_mat = dE2_mat + dE2_mat_cell;
end

% dE2_mat  = squeeze(sum(dE2_mat_cell,1));

end