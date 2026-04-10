function [hessian_mat, sorted_eigenvalues, sorted_eigenvectors] = get_qvm_hessian(V,C,KP,P0,KA,A0,strain)

% 2nd partial derivative of  energy with respect to vertex positions;
% \frac{\partial^2 E}{\partial h_{m\alpha} \partial h_{n,\beta}}
% see dE2_matrix.nb for analytical calculations

[Cell_Perim, Cell_Area] = get_cell_shape_LE(V,C,strain);
Ncell = length(C);
box_size = sqrt(Ncell);


Nvert = size(V,1);

hessian_mat = zeros(2*Nvert,2*Nvert);
for i = 1:Ncell % loop over cells
    hessian_cell = zeros(2*Nvert,2*Nvert);
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


	for m = 2:z+1
		for n = 2:z+1
			if mod(n,z) == mod(m,z)	
				%	m == n
				% x,x
				hessian_cell(neibs(m),neibs(n)) = ...
					2*KP(i)*(-P0(i) + Cell_Perim(i))*(-((Vx(m) - Vx(m+1))^2/((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2)^(3/2)) + 1/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) - (Vx(m) - Vx(m-1))^2/((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)^(3/2) + 1/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + ...
	 				2*KP(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2))^2 + (KA(i)*(Vy(m+1) - Vy(m-1))^2)/2;

				% y, y
				hessian_cell(neibs(m)+Nvert,neibs(n)+Nvert) = ...
					(KA(i)*(Vx(m+1) - Vx(m-1))^2)/2 + 2*KP(i)*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2))^2 + ...
	 				2*KP(i)*(-P0(i) + Cell_Perim(i))*(1/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) - (Vy(m) - Vy(m+1))^2/((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2)^(3/2) + 1/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2) - (Vy(m) - Vy(m-1))^2/((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)^(3/2));

				% x, y == y,x
				hessian_cell(neibs(m),neibs(n)+Nvert) = ...
					2*KP(i)*(-P0(i) + Cell_Perim(i))*(-(((Vx(m) - Vx(m+1))*(Vy(m) - Vy(m+1)))/((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2)^(3/2)) - ((Vx(m) - Vx(m-1))*(Vy(m) - Vy(m-1)))/((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)^(3/2)) + ...
	 				2*KP(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2))*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (KA(i)*(-Vx(m+1) + Vx(m-1))*(Vy(m+1) - Vy(m-1)))/2;

				hessian_cell(neibs(m)+Nvert,neibs(n)) = hessian_cell(neibs(m),neibs(n)+Nvert);

			elseif mod(n,z) == mod(m+1,z) 
				% n == m+1
				% x,x
				hessian_cell(neibs(m),neibs(n)) = ...
					(-2*KP(i)*(-P0(i) + Cell_Perim(i))*(Vy(m) - Vy(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 2*KP(i)*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	  				((Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (KA(i)*(-Vy(m) + Vy(n+1))*(Vy(n) - Vy(m-1)))/2;


				% y, y
				hessian_cell(neibs(m)+Nvert,neibs(n)+Nvert) = ...
					(KA(i)*(Vx(m) - Vx(n+1))*(-Vx(n) + Vx(m-1)))/2 - (2*KP(i)*(-P0(i) + Cell_Perim(i))*(Vx(m) - Vx(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 2*KP(i)*((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	  				((Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2));

				% x, y 
				hessian_cell(neibs(m),neibs(n)+Nvert) = ...
					(-((2*A0(i) - 2*Cell_Area(i))*KA(i)) + (4*KP(i)*(-P0(i) + Cell_Perim(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 4*KP(i)*((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	   				((Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + KA(i)*(Vx(m) - Vx(n+1))*(Vy(n) - Vy(m-1)))/2;

				% y, x
				hessian_cell(neibs(m)+Nvert,neibs(n)) = ...
					((2*A0(i) - 2*Cell_Area(i))*KA(i) + (4*KP(i)*(-P0(i) + Cell_Perim(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + KA(i)*(Vx(n) - Vx(m-1))*(Vy(m) - Vy(n+1)) + 4*KP(i)*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))* ...
	   				((Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)))/2;


			elseif mod(n,z) == mod(m-1,z) 
				% x,x
				hessian_cell(neibs(m),neibs(n)) = ...
					2*KP(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2))*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2)) - ...
	 				(2*KP(i)*(-P0(i) + Cell_Perim(i))*(Vy(m) - Vy(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + (KA(i)*(Vy(m) - Vy(n-1))*(Vy(m+1) - Vy(n)))/2;

				% y, y
				hessian_cell(neibs(m)+Nvert,neibs(n)+Nvert) = ...
					(KA(i)*(Vx(m) - Vx(n-1))*(Vx(m+1) - Vx(n)))/2 - (2*KP(i)*(-P0(i) + Cell_Perim(i))*(Vx(m) - Vx(n))^2)/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + 2*KP(i)*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2))* ...
					((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2));

				% x, y 
				hessian_cell(neibs(m),neibs(n)+Nvert) = ...
					((2*A0(i) - 2*Cell_Area(i))*KA(i) + (4*KP(i)*(-P0(i) + Cell_Perim(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2) + KA(i)*(-Vx(m) + Vx(n-1))*(Vy(m+1) - Vy(n)) + 4*KP(i)*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2))* ...
	   				((-Vy(m) + Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2)))/2;

				% y, x
				hessian_cell(neibs(m)+Nvert,neibs(n)) = ...
					(-((2*A0(i) - 2*Cell_Area(i))*KA(i)) + KA(i)*(-Vx(m+1) + Vx(n))*(Vy(m) - Vy(n-1)) + 4*KP(i)*((-Vx(m) + Vx(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2) + (-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2))* ...
	   				((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(n))/sqrt((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)) + (4*KP(i)*(-P0(i) + Cell_Perim(i))*(Vx(m) - Vx(n))*(Vy(m) - Vy(n)))/((Vx(m) - Vx(n))^2 + (Vy(m) - Vy(n))^2)^(3/2))/2;
				

			else % m and n are not neighboring vertices 
				% x,x
				hessian_cell(neibs(m),neibs(n)) = ...
					2*KP(i)*((-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (KA(i)*(-Vy(n-1) + Vy(n+1))*(Vy(m+1) - Vy(m-1)))/2;


				% y, y
				hessian_cell(neibs(m)+Nvert,neibs(n)+Nvert) = ...
					(KA(i)*(Vx(n-1) - Vx(n+1))*(-Vx(m+1) + Vx(m-1)))/2 + 2*KP(i)*((-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2));


				% x, y 
				hessian_cell(neibs(m),neibs(n)+Nvert) = ...
					2*KP(i)*((-Vy(n-1) + Vy(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vy(n) - Vy(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vx(m) - Vx(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vx(m) - Vx(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2)) + (KA(i)*(Vx(n-1) - Vx(n+1))*(Vy(m+1) - Vy(m-1)))/2;


				% y, x
				hessian_cell(neibs(m)+Nvert,neibs(n)) = ...
					(KA(i)*(Vx(m+1) - Vx(m-1))*(Vy(n-1) - Vy(n+1)))/2 + 2*KP(i)*((-Vx(n-1) + Vx(n))/sqrt((Vx(n-1) - Vx(n))^2 + (Vy(n-1) - Vy(n))^2) + (Vx(n) - Vx(n+1))/sqrt((Vx(n) - Vx(n+1))^2 + (Vy(n) - Vy(n+1))^2))*((Vy(m) - Vy(m+1))/sqrt((Vx(m) - Vx(m+1))^2 + (Vy(m) - Vy(m+1))^2) + (Vy(m) - Vy(m-1))/sqrt((Vx(m) - Vx(m-1))^2 + (Vy(m) - Vy(m-1))^2));

				
			end
		end
    end


    hessian_mat = hessian_mat + hessian_cell;
end


diag_vals = diag(hessian_mat);
hessian_mat = hessian_mat - diag(diag(hessian_mat));
hessian_mat = 0.5*(hessian_mat+hessian_mat');
LL = size(hessian_mat,1);
hessian_mat(1:(LL+1):(LL)^2) = diag_vals;
[eigenvectors, eigenvalues] = eig(hessian_mat);
[sorted_eigenvalues,ind] = sort(diag(eigenvalues));
sorted_eigenvectors = eigenvectors(:, ind);

end