% Periodic Boundary Conditions version
function [perim_list, area_list] = get_cell_shape(V,newC)

        Npts = size(newC,1);
        box_size = sqrt(Npts);
        area_list = zeros(Npts,1);
        perim_list = zeros(Npts,1);
        for i=1:Npts
            X = V(newC{i},:);
            Z = size(X,1);
            Ximg = zeros(Z,2);
            Ximg(1,:) = X(1,:);
            Ximg(2:Z,1) = X(1,1) + pbc_dist(X(2:Z,1)-X(1,1),box_size);
            Ximg(2:Z,2) = X(1,2) + pbc_dist(X(2:Z,2)-X(1,2),box_size);

%             [~, area_list(i)] = convhull(Ximg(1:Z,1),Ximg(1:Z,2)); 
%             
%             dx = diff([Ximg(:,1); Ximg(1,1)]);
%             dy = diff([Ximg(:,2); Ximg(1,2)]);
%             perim_list(i)=sum(sqrt(dx.^2+dy.^2));
            
            
            dx = Ximg( [ 2:Z 1 ], 1) - Ximg(:,1);
            dy = Ximg( [ 2:Z 1 ], 2) - Ximg(:,2);

            % summations for CW boundary integrals
            area_list(i) = abs(sum( Ximg(:,2).*dx - Ximg(:,1).*dy )/2);
            
            perim_list(i) = sum( sqrt( dx.*dx +dy.*dy ) );
        end
%         if any(area_list<0)
%             disp('some areas are negative!!!');
%         end
end
