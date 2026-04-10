function [vvx,vvy] = get_vertex_velo_viscous(V,topo,VF,B,strain)
% %%%%% VF is active force + conserved force, B is viscous coefficient
    N = size(V,1)/2;
    coeff_mat = zeros(4*N,4*N);
    box_size = sqrt(N);
    adjVV = topo.VVadj;
    Vx =V(:,1);
    Vy =V(:,2);
    neib_Vx = adjVV.*Vx';
    neib_Vy = adjVV.*Vy';
    xi_xj = -(neib_Vx-Vx.*adjVV);              %%%% Vxi - Vxj matrix, where vertices j are neib to vertex i
    yi_yj = -(neib_Vy-Vy.*adjVV);              %%%% Vyi - Vyj matrix, where vertices j are neib to vertex i
    [ex,ey] = pbc_dist_LE(xi_xj,yi_yj,box_size,strain);
    ex=-ex;
    ey=-ey;

    edge_L  = sqrt(ex.^2+ey.^2)   ;        %%%%% edge length matrix correspond to adjVV
    ex_hat = ex./edge_L;
    ex_hat(isnan(ex_hat))=0;             %%%%% x-comp unit vector of hi - hj correspond to adjvv
    ey_hat = ey./edge_L;
    ey_hat(isnan(ey_hat))=0;             %%%%% y-comp unit vector of hi - hj correspond to adjvv
    coef_xx = eye(size(adjVV)).*(B*sum(ex.*ex_hat,2)+1)-B*ex.*ex_hat;
    coef_xy =  eye(size(adjVV)).*(B*sum(ey.*ex_hat,2))-B*ey.*ex_hat;
    coef_yx = eye(size(adjVV)).*(B*sum(ey.*ex_hat,2))-B*ex.*ey_hat;
    coef_yy = eye(size(adjVV)).*(B*sum(ey.*ey_hat,2)+1)-B*ey.*ey_hat;

    %%%%%%%%%%%%%%%%%% fill the coef matrix
    coeff_mat(1:2*N,1:2*N)=coef_xx;
    coeff_mat(1:2*N,2*N+1:4*N)=coef_xy;
    coeff_mat(2*N+1:4*N,1:2*N)=coef_yx;
    coeff_mat(2*N+1:4*N,2*N+1:4*N)=coef_yy;
    A = [VF(:,1) ; VF(:,2)];
    vv = (coeff_mat^-1)*A;
    vvx = vv(1:2*N);
    vvy = vv(2*N+1:end);  
end