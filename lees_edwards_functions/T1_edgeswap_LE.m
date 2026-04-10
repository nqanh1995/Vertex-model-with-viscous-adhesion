function [V,C,topo] = T1_edgeswap_LE(V,C,topo,T1_edge_ID,strain)

box_size = sqrt(topo.Ncell);

u0 = topo.E_list(T1_edge_ID,1);
v0 = topo.E_list(T1_edge_ID,2);


temp = find(topo.CEadj(:,T1_edge_ID));
c1 = temp(1);
c2 = temp(2);
zz = length(C{c1});
ind1 = find(C{c1} == u0);
if C{c1}((1 + mod(ind1+1-1, zz))) == v0
else
    c1 = temp(2);
    c2 = temp(1);
end


z1 = length(C{c1});
z2 = length(C{c2});

if z1 < 4 | z2 < 4
    disp('there are triangles, no T1 possible!');
    return;
end


cell_neibs = topo.CVadj(:,u0);
cell_neibs([c1,c2]) = 0;
cu = find(cell_neibs);

cell_neibs = topo.CVadj(:,v0);
cell_neibs([c1,c2]) = 0;
cv = find(cell_neibs);



cell_list_cu = C{cu};
zz = length(cell_list_cu);
ind0 = find(cell_list_cu == u0);
ind1 = (1 + mod(ind0-1-1, zz));
ind2 = (1 + mod(ind0+1-1, zz));
u1 = cell_list_cu(ind1);
u2 = cell_list_cu(ind2);
% T1 update cell list for cell cu
cell_list_cu = [cell_list_cu(1:ind0),v0,cell_list_cu(ind0+1:end)];


cell_list_cv = C{cv};
zz = length(cell_list_cv);
ind0 = find(cell_list_cv == v0);
ind1 = (1 + mod(ind0-1-1, zz));
ind2 = (1 + mod(ind0+1-1, zz));
v1 = cell_list_cv(ind1);
v2 = cell_list_cv(ind2);
% T1 update cell list for cell cu
cell_list_cv = [cell_list_cv(1:ind0),u0,cell_list_cv(ind0+1:end)];

cell_list_c1 = C{c1};
% T1 update cell list for cell c1
cell_list_c1(cell_list_c1 == u0) = [];


cell_list_c2 = C{c2};
% T1 update cell list for cell c2
cell_list_c2(cell_list_c2 == v0) = [];


%% update the data structure
% cell list
C{cu} = cell_list_cu;
C{cv} = cell_list_cv;
C{c1} = cell_list_c1;
C{c2} = cell_list_c2;
% % v-v adj
topo.VVadj(u0,u2) = 0; topo.VVadj(u2,u0) = 0;
topo.VVadj(v0,v2) = 0; topo.VVadj(v2,v0) = 0;
topo.VVadj(u2,v0) = 1; topo.VVadj(v0,u2) = 1;
topo.VVadj(u0,v2) = 1; topo.VVadj(v2,u0) = 1;
% % c-v adj
topo.CVadj(c1,u0) = 0;
topo.CVadj(c2,v0) = 0;
topo.CVadj(cu,v0) = 1;
topo.CVadj(cv,u0) = 1;
% % edges
E1 = find(...
    (topo.E_list(:,1) == u0 & topo.E_list(:,2) == u2) | ...
    (topo.E_list(:,2) == u0 & topo.E_list(:,1) == u2));
E4 = find(...
    (topo.E_list(:,1) == v0 & topo.E_list(:,2) == v2) | ...
    (topo.E_list(:,2) == v0 & topo.E_list(:,1) == v2));
% % updating edge list
topo.E_list(E1,:) = [v0,u2];
topo.E_list(E4,:) = [u0,v2];
% % v-e adj
topo.VEadj(u0,E1) = 0; topo.VEadj(v0,E1) = 1;
topo.VEadj(v0,E4) = 0; topo.VEadj(u0,E4) = 1;
% % c-e adj
topo.CEadj(c1,T1_edge_ID) = 0;
topo.CEadj(c2,T1_edge_ID) = 0;
topo.CEadj(cu,T1_edge_ID) = 1;
topo.CEadj(cv,T1_edge_ID) = 1;
% rotate the vertices of the T1 edge

ex = V(v0,1) - V(u0,1);
ey = V(v0,2) - V(u0,2);
[ex,ey] = pbc_dist_LE(ex,ey,box_size,strain);


xv0 = V(u0,1) + ex;
yv0 = V(u0,2) + ey;
xu0 = V(u0,1);
yu0 = V(u0,2);
xmid = 0.5*(xv0 +  xu0);
ymid = 0.5*(yv0 +  yu0);


V(v0,1) = -(yv0-ymid)+xmid;
V(v0,2) = (xv0 - xmid)+ymid;
V(u0,1) = -(yu0-ymid)+xmid;
V(u0,2) = (xu0 - xmid)+ymid;
V = position_mod_LE(V,box_size,strain);


% V(v0,1) = mod(-(yv0-ymid)+xmid,box_size);
% V(v0,2) = mod((xv0 - xmid)+ymid,box_size);
% V(u0,1) = mod(-(yu0-ymid)+xmid,box_size);
% V(u0,2) = mod((xu0 - xmid)+ymid,box_size);



end