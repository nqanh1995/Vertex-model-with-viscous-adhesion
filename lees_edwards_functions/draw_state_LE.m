function ff = draw_state_LE(V,C,strain,fignum,c,lw,opac)

if fignum == 0

    ff = figure();
else
    ff = figure(fignum);
end
N = size(C,1);
box_size = sqrt(N);

cm = colormap(jet(N));
hold on;
for i=1:N
    X = V(C{i},1);
    Y = V(C{i},2);
    z = size(X,1);
    % X = [X(:);X(1)];
    % Y = [Y(:);Y(1)];

    % for j = 1:z
    %     dx = X(j+1)-X(j);
    %     dy = Y(j+1)-Y(j);
    %     [dx,dy] = pbc_dist_LE(dx,dy,box_size,strain);
    % 
    % 
    % 
    %     plot([X(j),X(j)+dx],[Y(j),Y(j)+dy],'-k');
    % 
    % end
    dx = X(2:z)-X(1);
    dy = Y(2:z)-Y(1);
    [dx,dy] = pbc_dist_LE(dx,dy,box_size,strain);

    X(2:z) = X(1) + dx;
    Y(2:z) = Y(1) + dy;


    h = plot([X(:);X(1)],[Y(:);Y(1)],color=c,LineWidth=lw);
    alpha(opac);
    % h = patch([X(:);X(1)],[Y(:);Y(1)],'w','facecolor',cm(i,:));

    hold on;
end
axis([-2,box_size+2,-2,box_size+2])
rectangle('Position',[0,0,box_size,box_size])
pbaspect([1,1,1]);
set(gca, 'Color', 'none'); % Sets axes background
% text(box_size/2,-1,['strain = ',num2str(strain)])

end