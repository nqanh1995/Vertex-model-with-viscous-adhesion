function ff = draw_state(V,C)
ff = figure();
N = size(C,1);
box_size = sqrt(N);


hold on;
for i=1:N
    X = V(C{i},:);
    Z = size(X,1);
    Ximg = zeros(Z,2);
    Ximg(1,:) = X(1,:);

    Ximg(2:Z,1) = X(1,1) + pbc_dist(X(2:Z,1)-X(1,1),box_size);
    Ximg(2:Z,2) = X(1,2) + pbc_dist(X(2:Z,2)-X(1,2),box_size);
    meanxy = mean(Ximg,1);
    if meanxy(1) - X(1,1) > 0.5*box_size
        Ximg(:,1) = Ximg(:,1) - box_size;
    end
    if meanxy(1) - X(1,1) < -0.5*box_size
        Ximg(:,1) = Ximg(:,1) + box_size;
    end
    if meanxy(2) - X(1,2) > 0.5*box_size
        Ximg(:,2) = Ximg(:,2) - box_size;
    end
    if meanxy(2) - X(1,2) < -0.5*box_size
        Ximg(:,2) = Ximg(:,2) + box_size;
    end



    h = plot([Ximg(:,1);Ximg(1,1)],[Ximg(:,2);Ximg(1,2)],'-k');
    % if Z == 6
    %     h = patch([Ximg(:,1);Ximg(1,1)],[Ximg(:,2);Ximg(1,2)],'w');
    % elseif Z == 7
    %     h = patch([Ximg(:,1);Ximg(1,1)],[Ximg(:,2);Ximg(1,2)],'c');
    % elseif Z == 5
    %     h = patch([Ximg(:,1);Ximg(1,1)],[Ximg(:,2);Ximg(1,2)],'y');
    % else 
    %     h = patch([Ximg(:,1);Ximg(1,1)],[Ximg(:,2);Ximg(1,2)],'k');
    % 
    % end

    hold on;
end
axis([-2,box_size+2,-2,box_size+2])
rectangle('Position',[0,0,box_size,box_size])
pbaspect([1,1,1]);
set(gca, 'Color', 'none'); % Sets axes background

end