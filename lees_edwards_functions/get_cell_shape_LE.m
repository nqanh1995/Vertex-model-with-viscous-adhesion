% Periodic Boundary Conditions version
function [plist, alist] = get_cell_shape_LE(V,C,strain)

        N = size(C,1);
        box_size = sqrt(N);
        alist = zeros(N,1);
        plist = zeros(N,1);

        for i=1:N
            X = V(C{i},1);
            Y = V(C{i},2);
            z = size(X,1);

            dx = X(2:z)-X(1);
            dy = Y(2:z)-Y(1);
            [dx,dy] = pbc_dist_LE(dx,dy,box_size,strain);

            X(2:z) = X(1) + dx;
            Y(2:z) = Y(1) + dy;

            dx = X( [ 2:z, 1 ]) - X;
            dy = Y( [ 2:z, 1 ]) - Y;

            alist(i) = abs(sum( Y.*dx - X.*dy )/2);
            
            plist(i) = sum( sqrt( dx.*dx +dy.*dy ) );

        end


end
