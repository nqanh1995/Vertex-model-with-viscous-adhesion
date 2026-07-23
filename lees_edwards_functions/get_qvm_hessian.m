function [hessian_mat, sorted_eigenvalues, sorted_eigenvectors] = get_qvm_hessian(V,C,KP,P0,KA,A0,strain)
% GET_QVM_HESSIAN  Build the vertex-model energy Hessian and return its
% eigen-decomposition, under Lees-Edwards (sheared) periodic boundary
% conditions.
%
%   [hessian_mat, sorted_eigenvalues, sorted_eigenvectors] =
%   GET_QVM_HESSIAN(V, C, KP, P0, KA, A0, strain) assembles the full
%   system Hessian d^2E/(dh_m dh_n) exactly as in get_dE2_mat_vm_LE.m
%   (see that file for a detailed explanation of the four per-vertex-pair
%   cases below), then symmetrizes it and computes its eigenvalues and
%   eigenvectors, sorted in ascending order. The eigenvalues characterize
%   the tissue's linear response modes (e.g. soft/zero modes indicate
%   marginal mechanical stability).
%
% Inputs:
%   V, C     - vertex positions and cell-vertex lists
%   KP, P0   - perimeter elasticity coefficient and target perimeter, per
%              cell
%   KA, A0   - area elasticity coefficient and target area, per cell
%   strain   - current Lees-Edwards shear strain
%
% Outputs:
%   hessian_mat         - 2*Nvert x 2*Nvert symmetrized Hessian matrix
%   sorted_eigenvalues  - eigenvalues of hessian_mat, ascending order
%   sorted_eigenvectors - corresponding eigenvectors (columns), in the
%                         same order as sorted_eigenvalues
%
% Note: the per-vertex-pair closed-form expressions below are derived
% symbolically (see dE2_matrix.nb) and reproduced as-is; do not hand-edit.

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
    % Unwrap this cell's vertices relative to its last vertex, so all
    % vertices are expressed in the same (sheared) periodic image before
    % computing the derivatives below.
    [dx,dy] = pbc_dist_LE(Vx-Vx(end),Vy-Vy(end),box_size,strain);
    Vx = Vx(end) + dx;
    Vy = Vy(end) + dy;

    % Extend the vertex list cyclically so each vertex has a well-defined
    % previous/next neighbor around the polygon, including wrap-around.
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];

    % Loop over every pair of vertices (m, n) belonging to this cell; the
    % applicable case depends on whether n is the same vertex as m, the
    % next/previous vertex around the polygon, or unrelated to m.
    for m = 2:z+1
        for n = 2:z+1
            if mod(n,z) == mod(m,z) 
                % Case 1: n == m
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
                % Case 2: n is the vertex immediately after m
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
                % Case 3: n is the vertex immediately before m
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
                

            else
                % Case 4: m and n are not adjacent (not connected by an edge)
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

    % Add this cell's Hessian contribution into the full system matrix.
    hessian_mat = hessian_mat + hessian_cell;
end

% Symmetrize the Hessian: the analytical off-diagonal blocks (x,y) and
% (y,x) for the "n==m" case are set equal by construction, but small
% floating-point asymmetries can remain elsewhere, so the diagonal is
% preserved exactly and the off-diagonal part is explicitly symmetrized.
diag_vals = diag(hessian_mat);
hessian_mat = hessian_mat - diag(diag(hessian_mat));
hessian_mat = 0.5*(hessian_mat+hessian_mat');
LL = size(hessian_mat,1);
hessian_mat(1:(LL+1):(LL)^2) = diag_vals;   % restore the exact diagonal

% Eigen-decompose the (now symmetric) Hessian and sort ascending, so
% sorted_eigenvalues(1) is the softest mode of the tissue.
[eigenvectors, eigenvalues] = eig(hessian_mat);
[sorted_eigenvalues,ind] = sort(diag(eigenvalues));
sorted_eigenvectors = eigenvectors(:, ind);

end