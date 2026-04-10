function ff = draw_edges_LE(V,topo,strain,fignum)

if fignum == 0

    ff = figure();
else
    ff = figure(fignum);
end

box_size = sqrt(topo.Ncell);

cm = colormap(hsv(topo.Nedge));
hold on;
for i = 1:topo.Nedge


    v1 = topo.E_list(i,1);
    v2 = topo.E_list(i,2);
    dx = V(v2,1)-V(v1,1);
    dy = V(v2,2)-V(v1,2);
    [dx,dy] = pbc_dist_LE(dx,dy,box_size,strain);

    plot([V(v1,1),V(v1,1)+dx],[V(v1,2),V(v1,2)+dy],'-','color',cm(i,:),'LineWidth',2);

end
axis([-2,box_size+2,-2,box_size+2])
rectangle('Position',[0,0,box_size,box_size])
pbaspect([1,1,1]);
set(gca, 'Color', 'none'); % Sets axes background
text(box_size/2,-1,['strain = ',num2str(strain)])

end