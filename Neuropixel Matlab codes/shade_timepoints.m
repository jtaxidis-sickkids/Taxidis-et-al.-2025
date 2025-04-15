function shade_timepoints(lims,tp)

Y = lims([1 2 2 1]);
patch([tp(1)*[1 1] tp(2)*[1 1]], Y,[.5 .5 .5],'EdgeColor','none', 'FaceAlpha',0.2);
if length(tp) > 2
    patch([tp(3)*[1 1], tp(4)*[1 1]], Y,[.5 .5 .5],'EdgeColor','none', 'FaceAlpha',0.2);
    patch([tp(5)*[1 1], tp(6)*[1 1]], Y,[.5 .5 .5],'EdgeColor','none', 'FaceAlpha',0.2);
end