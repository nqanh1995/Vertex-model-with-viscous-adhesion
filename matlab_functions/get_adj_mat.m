function topo =  get_adj_mat(V,C)
%%%%%%%%%% C is cell list
%%%%%%%%%% V is vertex position
Ncell = size(C,1);
Nvert = size(V,1);

% obtain the vertex-vertex adjacency matrix
VVadj = zeros(Nvert,Nvert);
for i = 1:Ncell
    idx = C{i};
    ZC = length(idx);
    for j = 1:(ZC-1)
        VVadj(idx(j),idx(j+1)) = 1;
    end
    VVadj(idx(ZC),idx(1)) = 1;
end
VVadj = VVadj + transpose(VVadj);
VVadj = logical(VVadj);

%  edge list

[vv1,vv2] = find(triu(VVadj));
E_list = [vv1,vv2];
Nedge = size(E_list,1);

% obtain the vertex-edge adjacency matrix
VEadj = zeros(Nvert,Nedge); 

for eid = 1:Nedge
    VEadj(vv1(eid),eid) = 1;
    VEadj(vv2(eid),eid) = 1;
    
end


%  Cell-Vertex adjacency matrix 
CVadj = zeros(Ncell,Nvert);
for i = 1:Ncell
    CVadj(i,C{i}) = 1;
end



% obtain the  Cell-edge adjacency matrix
CEadj = zeros(Ncell,Nedge);
for i = 1:Ncell
    zz = length(C{i});
    for j = 1:(zz-1)
        eid = intersect(find(VEadj(C{i}(j),:)),find(VEadj(C{i}(j+1),:)));
        CEadj(i,eid) = 1;
        
    end
    eid = intersect(find(VEadj(C{i}(zz),:)),find(VEadj(C{i}(1),:)));
    CEadj(i,eid) = 1;
end
CEadj = logical(CEadj);


topo =  struct(...
                  'Ncell',Ncell,...
                  'Nvert',Nvert,...
                  'Nedge',Nedge,...
                  'VVadj',VVadj,...
                  'CEadj',CEadj,...
                  'CVadj',CVadj,...
                  'E_list',E_list,...
                  'VEadj',VEadj...
                  );


end
