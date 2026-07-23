function dE2_mat = get_dE2_mat_vm_LE(V,newC,K_A,K_P,A0,P0,strain)
% GET_DE2_MAT_VM_LE  Assemble the full Hessian (second derivative) matrix
% of the vertex-model energy with respect to vertex positions, under
% Lees-Edwards (sheared) periodic boundary conditions.
%
%   dE2_mat = GET_DE2_MAT_VM_LE(V, newC, K_A, K_P, A0, P0, strain)
%   computes, for every cell, the analytical second derivatives
%       d^2 E / (d h_{m,alpha} d h_{n,beta})
%   of that cell's area+perimeter elastic energy with respect to pairs of
%   its own vertices (m, n), for both coordinate components (alpha, beta
%   in {x, y}), and sums these per-cell contributions into the full
%   2*Nvert x 2*Nvert system Hessian. This is used for linear stability /
%   shear-modulus calculations (see get_shear_modulus_LE.m,
%   get_qvm_hessian.m).
%
% The closed-form expressions below were derived symbolically (see the
% companion Mathematica notebook dE2_matrix.nb) and are reproduced here
% as-is; they are algebraically dense but not conceptually complex --
% each block below just fills in one of four cases, depending on whether
% vertex n is: the same vertex as m, the next vertex around the cell
% polygon, the previous vertex, or unrelated (not adjacent) to m. Do not
% hand-edit these expressions; regenerate from the notebook if the energy
% functional changes.
%
% Inputs:
%   V      - Nvert x 2 array of vertex positions
%   newC   - Ncell x 1 cell array of vertex indices (in polygon order)
%            belonging to each cell
%   K_A,A0 - area elasticity coefficient and target area, per cell
%   K_P,P0 - perimeter elasticity coefficient and target perimeter, per
%            cell
%   strain - current Lees-Edwards shear strain
%
% Output:
%   dE2_mat - 2*Nvert x 2*Nvert Hessian matrix; rows/columns 1:Nvert are
%             the x-components, Nvert+1:2*Nvert are the y-components

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

    % Unwrap this cell's vertices relative to its first vertex, so all
    % vertices are expressed in the same (sheared) periodic image before
    % computing the finite-difference-style derivatives below.
    dx = Vx(2:z)-Vx(1);
    dy = Vy(2:z)-Vy(1);
    [dx,dy]=pbc_dist_LE(dx,dy,box_size,strain);
    Vx(2:z) = Vx(1) + dx;
    Vy(2:z) = Vy(1) + dy;

    % Extend the vertex list cyclically so each vertex has a well-defined
    % previous/next neighbor around the polygon, including wrap-around.
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];

    % Loop over every pair of vertices (m, n) belonging to this cell, and
    % fill in the corresponding Hessian block. Which closed-form
    % expression applies depends on the topological relationship between
    % m and n around the polygon (same vertex, adjacent, or unrelated).
    for m = 2:z+1
        for n = 2:z+1
            if mod(n,z) == mod(m,z) 
                % Case 1: n == m (diagonal block -- second derivative of
                % the cell's energy with respect to the same vertex twice).
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
                % Case 2: n is the vertex immediately after m around the
                % polygon (m and n share an edge).
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
                % Case 3: n is the vertex immediately before m around the
                % polygon (m and n share an edge, opposite orientation of
                % Case 2).
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
                

            else
                % Case 4: m and n are not adjacent (no shared edge) -- the
                % only coupling between them is through the cell's area
                % term (perimeter terms only couple adjacent vertices).
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

    % Add this cell's Hessian contribution into the full system matrix.
    dE2_mat = dE2_mat + dE2_mat_cell;
end
end