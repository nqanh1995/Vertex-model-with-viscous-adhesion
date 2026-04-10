% Periodic Boundary Conditions version
function [perimeters, areas, a1, a2, orientation] = cell_shape_info(V,newC)

        Npts = size(newC,1);
        box_size = sqrt(Npts);
        areas = zeros(Npts,1);
        perimeters = zeros(Npts,1);
        orientation = zeros(Npts,1);
        a1 = zeros(Npts,1);
        a2 = zeros(Npts,1);
        for i=1:Npts
            vx = V(newC{i},1);
            vy = V(newC{i},2);
            z = size(vx,1);
            vx(2:z) = vx(1) + pbc_dist(vx(2:z)-vx(1),box_size);
            vy(2:z) = vy(1) + pbc_dist(vy(2:z)-vy(1),box_size);
            ccw_order=vertices_to_cw_polygon(vx,vy);
            vx = vx(ccw_order);
            vy = vy(ccw_order);


            [geom,~,cpmo] = polygeom(vx,vy);
            areas(i) = geom(1);
            perimeters(i) = geom(4);
            a1(i) = cpmo(3);
            a2(i) = cpmo(1);
            orientation(i) = cpmo(2);

        end
    
end
