function scatter_plot_comparison(X,Y,p)

hold on;

if isrow(X)
    X = X';
end
if isrow(Y)
    Y = Y';
end

scatter(X,Y,8,'k');%,'filled');

line([min([X;Y]),max([X;Y])],[min([X;Y]),max([X;Y])],'Color','k');

mx = nanmean(X);
my = nanmean(Y);
sx = SEM(X);
sy = SEM(Y);

errorbar(mx,my,sy,sy,sx,sx,'sk','Markerfacecolor','k','Markersize',10);

if p < 0.05
    plot(mx, my*1.2,'*k','MarkerSize',6);
end

axis tight