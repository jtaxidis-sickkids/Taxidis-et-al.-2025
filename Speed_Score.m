function Speed_Score(INclass)

%% LOAD POOLED DATA
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'R','MOT','bins','timepoints');

%% STACK IN SINGLE COLUMNS
R = vertcat(R{:});
MOT = vertcat(MOT{:});

ls = length(R); % Number of sessions

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
Rall = [];
MOTall = {};

for i = 1:ls
    for r = 1:size(R{i},1)
        Rall{count,1} = squeeze(R{i}(r,:,:));                               % Keep only that cell's rates [time x trials]
        MOTall{count,1} = MOT{i};                                           % (Same dimensions)
        count = count + 1;
    end
end
R = Rall;
MOT = MOTall;
clear Rall MOTall

ls = count-1;

%% SMOOTH RATES AND DOWNSAMPLE MOTION
for i = 1:ls
    MOT{i} = MOT{i}(1:100:end,:);                                           % Downsample to have same vector length as rates
    for tr = 1:size(R{i},2)
        R{i}(:,tr) = smooth(R{i}(:,tr),5,'moving');
    end
end

%% MAKE TRIAL VECTOR
lb = length(bins);
trial = zeros(lb,1);
od1 = (bins >= timepoints(1) & bins <= timepoints(2));
od2 = (bins >= timepoints(3) & bins <= timepoints(4));
ods = logical(od1 + od2);
trial(ods) = 1;                                                             % Vector of zeros with ones only during odors

%% GET LOCOMOTION-RATES CORRELATIONS (SPEED SCORE)
Corr = nan(ls,1);
for i = 1:ls
    C = corr(R{i},MOT{i},'type','Pearson');                                 % Correlation of rate and motion over concatenated trials
    C = diag(C);                                                            % Keep only correlations of each trial with itself
    Corr(i) = mean(C);
end

%% GET TRIAL-RATES CORRELATIONS (ODOR SCORE)
CorrOD = nan(ls,1);
for i = 1:ls
    C = corr(R{i},trial,'type','Pearson');                                  % Correlation of theta power and the trial-vector [time x 1]
    CorrOD(i) = mean(C);
end

%% REPEAT FOR SHUFFLES
reps = 1000;
Corr_sh = zeros(ls,reps);
CorrOD_sh = zeros(ls,reps);

rng(1);

for i = 1:ls
    Ntr = size(R{i},2);                                                     % Number of trials
    lags = 2*rand(reps,Ntr) - 1;                                            % 1000 random numbers in [-1 1] range
    lags = floor(lags*lb/2);                                                % Get random (round) lags up to +- duration/2 of total time
    for r = 1:reps
        Rsh = nan(lb,Ntr);
        for tr = 1:Ntr
            Rsh(:,tr) = circshift(R{i}(:,tr),lags(r,tr),1);                 % Shift time axis of each trial separately
        end
        C = corr(Rsh,MOT{i},'type','Pearson');
        C = diag(C);
        Corr_sh(i,r) = mean(C);
        
        C = corr(Rsh,trial,'type','Pearson');
        CorrOD_sh(i,r) = mean(C);
    end
end

%% FIND CELLS WITH SIGNIFICANT SPEED SCORES
maxthr = prctile(Corr_sh,95,2);                                             % 95% percentile for each cell over all repetitions
speedcells = Corr > maxthr;
speedcells = sum(speedcells)/ls*100;                                        % Ratio of significant speed cells

%% FIND CELLS WITH SIGNIFICANT ODOR SCORES
maxthr = prctile(CorrOD_sh,95,2);                                           % 95% percentile for each cell over all repetitions
odorcells = CorrOD > maxthr;
odorcells = sum(odorcells)/ls*100;                                          % Ratio of significant speed cells

%% PLOT
figure;
subplot(131); hold on
plot([-0.3:0.01:0.29], histcounts(Corr_sh(:),[-0.3:0.01:0.3])/reps,'r');    % Histogram of shuffles (normalized per shuffle)
histogram(Corr,[-0.3:0.01:0.3]);                                            % Histogram of normal correlations
% set(gca,'YScale','log')
title([num2str(speedcells),'% significant speed cells']);

%% PLOT
subplot(132); hold on
plot([-0.3:0.01:0.29], histcounts(CorrOD_sh(:),[-0.3:0.01:0.3])/reps,'r');
histogram(CorrOD,[-0.3:0.01:0.3]);
% set(gca,'YScale','log')
title([num2str(odorcells),'% significant odor cells']);

%% PLOT SPEED SCORE-ODOR SCORE
[~,p] = ttest(Corr,CorrOD);           % Paired sample ttest

subplot(133);
scatter_plot_comparison(Corr,CorrOD,p);
xlabel('Speed scores');
ylabel('Odor scores');
title(num2str(p));