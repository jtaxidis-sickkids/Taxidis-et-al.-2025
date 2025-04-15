function plot_mean_SE(data,pvalue,cols)

hold on;
for i = 1:length(data)
    M = nanmean(data{i});
    S = SEM(data{i});
    
    bar(i,M,'EdgeColor',cols(i,:),'FaceColor','w','LineWidth',2);
%     plot(i+randn(size(data{i}))*0.1, data{i},'o','Color',cols(i,:),'Markersize',3);
    errorbar(i,M,S,'k','Linewidth',2.5);
end

if ~isempty(pvalue)
    for i = 1:length(data)-1
        for j = i+1 : length(data)
            plot_significance(pvalue(i,j),i,j,max([nanmean(data{i}) nanmean(data{j})]),nanmean(data{i})*0.05);
        end
    end
end



%% -------------------------
function mb = plot_significance(p,loc1,loc2,mbold,incr)

if p <= 0.05
    mb = mbold + incr*4;
    ml = mean([loc1,loc2]);
    
    line([loc1,loc1] , [mb, mb+incr],'Color','k');
    line([loc2,loc2] , [mb, mb+incr],'Color','k');
    line([loc1,loc2] , [mb+incr, mb+incr],'Color','k');
    
    if p <= 0.05 && p > 0.01
        plot(ml , mb+incr*2,'*k','MarkerSize',5);
    elseif p <= 0.01 && p > 0.001
        plot(ml + [-ml/50,ml/50] , [mb+incr*2, mb+incr*2] , '*k','MarkerSize',5);
    elseif p <= 0.001
        plot(ml + [-ml/25,0,ml/25] , [mb+incr*2, mb+incr*2, mb+incr*2] , '*k','MarkerSize',5);
    end
else
    mb = mbold;
end

