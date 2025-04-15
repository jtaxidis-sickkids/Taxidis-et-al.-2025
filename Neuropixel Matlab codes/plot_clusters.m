function plot_clusters(X,WF,C1,C2,time,C,sub,tit)
subplot(2,4,sub); hold on;
plot(X(C1,1),X(C1,2),'o','Color','k','Markerfacecolor','k');                % Plot the large cluster first
plot(X(C2,1),X(C2,2),'o','Color','r','Markerfacecolor','r');
if ~isempty(C)
    plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',2);
end
xlabel('Peak-Trough distance');
ylabel('Peak-Trough ratio');
title(tit);

subplot(2,4,sub+4); hold on;
fill_plot(time,nanmean(WF(C1,:)),SEM(WF(C1,:)),'k')
plot(time,nanmean(WF(C1,:)),'k','linewidth',2);
fill_plot(time,nanmean(WF(C2,:)),SEM(WF(C2,:)),'r')
plot(time,nanmean(WF(C2,:)),'r','linewidth',2);
xlabel('time (ms)');
ylabel('Average waveforms');
lc1 = length(C1);
lc2 = length(C2);
title(['PYs = ',num2str(lc1),' (',num2str(100*lc1/(lc1+lc2)),'%). INs = ',num2str(lc2),' (',num2str(100*lc2/(lc1+lc2)),'%']);
axis tight
