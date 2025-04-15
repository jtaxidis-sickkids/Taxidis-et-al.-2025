function Figure_Rates_stats(INclass)  

%% LOAD POOLED DATA AND STACK
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'R','MC','timepoints','bins');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    R{1,5} = [];     % DFF{2,3} = [];
    MC{1,5} = [];      % MC{2,3} = [];
end
R = vertcat(R{:});
MC = vertcat(MC{:});

ls = length(R);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
Rall = {};
MCall = [];

for i = 1:ls
    for r = 1:size(R{i},1)
        Rall{count,1} = squeeze(R{i}(r,:,:))';                          % Keep only that cell's Rates [trials x time]
        MCall(count,:) = MC{i}(r,:);    
        count = count + 1;
    end
end
R = Rall;
MC = MCall;
clear DFFall MCall 

ls = count-1;

%% FIRING RATES AT TRIAL SEGMENTS
Rchange = nan(ls,4);   

for i = 1:ls                                                                % For each cell    
    baseR = R{i}(:,bins > 0.2 & bins < timepoints(1));                      % BASELINE rate
    
    od1R  = R{i}(:,bins >= timepoints(1) & bins <= timepoints(2));          % Odor 1 tates
    od2R  = R{i}(:,bins >= timepoints(3) & bins <= timepoints(4));          % Odor 2 rates
    deR   = R{i}(:,bins > timepoints(2) & bins < timepoints(3));            % Delay rates
    postR = R{i}(:,bins > timepoints(4) & bins < timepoints(6)-0.2);        % Post rates
    
    baseR = mean(baseR,2);                                                  % Average across bins
    od1R = mean(od1R,2);
    od2R = mean(od2R,2);
    deR = mean(deR,2);
    postR = mean(postR,2);
    
    baseR = mean(baseR);                                                    % Then average across trials
    od1R = mean(od1R);
    od2R = mean(od2R);
    deR = mean(deR);
    postR = mean(postR);
    
    Rchange(i,:) = 100*([od1R deR od2R postR] - baseR) / baseR;             % Compute relative change
end

%% KEEP ONLY MCELLS
% ismcell = ~isnan(MC(:,1));
% Rchange = Rchange(ismcell,:);

%% SIGNIFICANCE AGAINST BASELINE
pb = nan(4,1);
for i = 1:4
    [~,pb(i)] = ttest(Rchange(:,i),0,'alpha',0.05,'tail','right');
end
[~,pb] = fdr(pb);


%% SIGNIFICANCE IN DIFFERENCES
p = nan(3,1);
p(1) = ranksum(Rchange(:,1),Rchange(:,2));
p(2) = ranksum(Rchange(:,1),Rchange(:,3));
p(3) = ranksum(Rchange(:,1),Rchange(:,4));
[~,p] = fdr(p)


%% PLOT
figure;
A = table(Rchange(:,1),Rchange(:,2), Rchange(:,3),Rchange(:,4),'VariableNames',{'Odor1','Delay','Odor2','Post'});
violinplot(A);
hold on;

line([0,5],[0 0],'Color','r');
title(num2str(pb'));

plot_significance(p(1),1,2,Rchange(:,1),Rchange(:,2))
plot_significance(p(2),1,3,Rchange(:,1),Rchange(:,3))
plot_significance(p(3),1,4,Rchange(:,1),Rchange(:,4))

