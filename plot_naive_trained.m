function p_prepost = plot_naive_trained(D)

Dpre = D(:,1:2);
Dpost = D(:,3:end);

p_prepost = ranksum(Dpre(:),Dpost(:),'tail','both');

dp = table(Dpre,Dpost,'VariableNames',{'Naive','Trained'});
violinplot(dp);
hold on;
plot_significance(p_prepost,1,2,Dpre,Dpost);
title(num2str(p_prepost));
