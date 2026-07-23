function [sm,dE2_mat] = get_shear_modulus_LE(Npts,P0,V,newC,gamma,ka)
% GET_SHEAR_MODULUS_LE  Compute the athermal (zero-temperature) shear
% modulus of the tissue at a given strain, splitting it into affine and
% non-affine contributions, under Lees-Edwards (sheared) periodic
% boundary conditions.
%
%   [sm, dE2_mat] = GET_SHEAR_MODULUS_LE(Npts, P0, V, newC, gamma, ka)
%   computes:
%     - shear_mod_affine: the modulus assuming every vertex simply moves
%       affinely with the applied shear (no relaxation);
%     - shear_mod_nonaffine: the correction from vertices additionally
%       relaxing (non-affinely) in response to the shear, obtained by
%       solving the Hessian system dE2_mat \ dE2_gamma_V;
%   and returns their difference sm = affine - nonaffine, which is the
%   standard decomposition used to compute linear shear moduli in vertex
%   models (see e.g. the athermal quasi-static shear / AQS literature).
%
% Inputs:
%   Npts  - number of cells (also used to set the box size = sqrt(Npts))
%   P0    - target perimeter (applied uniformly to all cells)
%   V, newC - vertex positions and cell-vertex lists
%   gamma - shear strain at which to evaluate the modulus
%   ka    - area elasticity coefficient (applied uniformly to all cells)
%
% Outputs:
%   sm      - shear modulus (affine - non-affine contributions)
%   dE2_mat - the full energy Hessian at this configuration/strain (see
%             get_dE2_mat_vm_LE.m), returned for reuse/inspection

V_count = size(V,1); 
box_size = sqrt(Npts);
P0 = P0*ones(Npts,1);
K_A = ka*ones(Npts,1); % set modulus for area term to be zero
K_P = ones(Npts,1); % set modulus for area term to be one
A0 = ones(Npts,1); 

gamma = mod(gamma,1);

% Full energy Hessian (see get_dE2_mat_vm_LE.m) and current cell shapes
% at this configuration/strain.
dE2_mat = get_dE2_mat_vm_LE(V,newC,K_A,K_P,A0,P0,gamma);
[pe, ar] = get_cell_shape_LE(V,newC,gamma); %#ok<ASGLU> % ar kept for clarity, unused below
affine_part = zeros(Npts,1);
dE2_gamma_V = zeros(Npts,2*V_count);
zlist = zeros(Npts,1);

for i = 1:Npts % loop over cells
    neibs = newC{i};
    z = length(neibs);
    zlist(i) = z;
    Vx = V(neibs,1);
    Vy = V(neibs,2);

    % Unwrap this cell's vertices relative to its first vertex, so all
    % vertices are expressed in the same (sheared) periodic image before
    % computing the derivatives below.
    dx = Vx(2:z)-Vx(1);
    dy = Vy(2:z)-Vy(1);
    [dx,dy]=pbc_dist_LE(dx,dy,box_size,gamma);
    Vx(2:z) = Vx(1) + dx;
    Vy(2:z) = Vy(1) + dy;

    % Extend the vertex list cyclically so each vertex has a well-defined
    % previous/next neighbor around the polygon, including wrap-around.
    neibs = [neibs, neibs(1), neibs(2)];
    Vx = [Vx; Vx(1); Vx(2)];
    Vy = [Vy; Vy(1); Vy(2)];

    % Accumulate the pieces (part0, part1, part2) of the analytical
    % derivative of this cell's perimeter energy with respect to the
    % applied shear strain gamma, summed edge-by-edge around the polygon.
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

    % This cell's contribution to the affine shear modulus.
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
end

% Sum each cell's mixed-derivative contribution, then combine the affine
% modulus with the non-affine correction obtained by solving the Hessian
% system (dE2_mat \ dE2_gamma_V'): this is the standard projection of the
% affine-shear "forcing" onto the vertices' relaxation response.
dE2_gamma_V = sum(dE2_gamma_V,1);
shear_mod_affine = sum(affine_part)/Npts;
shear_mod_nonaffine = (dE2_gamma_V)*(dE2_mat\(dE2_gamma_V)')/Npts;
sm = shear_mod_affine - shear_mod_nonaffine;
end