function fill_plot(x,y,stdy,c)

hold on;
if size(x,1) > 1 && size(x,2) == 1                % If not given in row arrays
    x = x';
    y = y';    
    stdy = stdy';
end
fill([x, fliplr(x)],[y+stdy, fliplr(y-stdy)],c,'FaceAlpha',0.2,'EdgeColor','none');
plot(x,y,'-','Color',c,'Linewidth',1.5,'Markersize',4);