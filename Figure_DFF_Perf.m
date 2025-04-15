function Figure_DFF_Perf(INclass,cell_flag)

if nargin == 1, cell_flag = 'all'; end

%% LOAD POOLED DATA AND STACK IN SINGLE COLUMN
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'DFF','MC','PROG','time','timepoints');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    DFF{1,5} = [];    DFF{2,3} = [];
    MC{1,5} = [];   MC{2,3} = [];
    PROG{1,5} = [];   PROG{2,3} = [];
end

% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
if strcmp(INclass,'SOM')
    DFF(3,3:7) = DFF(3,4:8);
    MC(3,3:7) = MC(3,4:8);
    PROG(3,3:7) = PROG(3,4:8);
end

% KEEP ONLY TRAINED SESSIONS
DFF(:,1:2) = [];
MC(:,1:2) = [];
PROG(:,1:2) = [];

%% STACK DAYS
DFF = vertcat(DFF{:});
MC = vertcat(MC{:});
PROG = vertcat(PROG{:});

ls = length(DFF);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
DFFall = {};
MCall = [];
PROGall = {};

for i = 1:ls
    for r = 1:size(DFF{i},1)
        DFFall{count,1} = squeeze(DFF{i}(r,:,:))';                          % Keep only that cell's DFF [trials x time]
        MCall(count,:) = MC{i}(r,:);
        PROGall{count,1} = PROG{i}(:,2);
        count = count + 1;
    end
end
DFF = DFFall;
MC = MCall;
PROG = PROGall;
clear DFFall MCall PROGall 

%% CHOOSE CELL GROUP TO USE
Nr = size(MC,1);                            % Number of all cells
if strcmp(cell_flag,'all')
    k = true(Nr,1);
elseif strcmp(cell_flag,'mcells')
    k = ~isnan(MC(:,1));
elseif strcmp(cell_flag,'odorcells')
    k = MC(:,1) > 0;
elseif strcmp(cell_flag,'nonodorcells')
    k = MC(:,1) == 0;
elseif strcmp(cell_flag,'nonmcells')
    k = isnan(MC(:,1));
end

DFF = DFF(k);
PROG = PROG(k);

Nr = sum(k);

%% SMOOTH RAW DFF
for i = 1:Nr
    DFF{i} = sgolayfilt(DFF{i}, 1, 21, [] , 2);                             % Apply linear filter along time bins
end

%% GET HYPERPOLARIZATIONS   
minD = cell(Nr,1);
maxD = cell(Nr,1);
for c = 1:Nr                                                                % For each selected cell
    [minD{c},maxD{c}] = get_DFF_minmax(DFF{c},time,timepoints);             % NORMALIZE ODOR DFF TO BASELINE AND COMPUTE MIN/MAX
    minD{c} = abs(minD{c});                                                 % Keep absolute value of hyperpolarization
end

%% SPLIT CORRECT - ERROR TRIALS AND COMPUTE PERFORMANCE
for c = 1:Nr                                                                % For each selected cell
    PROG{c}(PROG{c} == 1 | PROG{c} == 4) = 1;                               % Set correct trials = 1
    PROG{c}(PROG{c} == 2 | PROG{c} == 3) = 0;                               % and errors = 0
end

Perf = 100*cellfun(@sum, PROG) ./ cellfun(@length, PROG);                   % Compute performance across recording of cell

%% COMPUTE HYPERPOL AMPLITUDE AND FREQUENCY OVER ALL TRIALS    
MminD = cellfun(@nanmean, minD);                                            % Get average hyperpolarization (non-significant ones are not accounted for)  
    
FminD = nan(Nr,1);
for c = 1:Nr                                                                % For each cell
    FminD(c) = 100 * sum(~isnan(minD{c})) / length(minD{c});                % Count ratio of significant hyperpolarizations over trials
end

%% COMPUTE HYPERPOL AMPLITUDE AND FREQUENCY OVER CORRECT VS ERROR
MminDc = nan(Nr,1);
MminDe = nan(Nr,1);
FminDc = nan(Nr,1);
FminDe = nan(Nr,1);

for c = 1:Nr                                                                % For each cell
    k = (PROG{c} == 1);                                                     % Find all correct trials of its video(s)
    MminDc(c) = nanmean(minD{c}(k));                                        % Get average hyperpolarization for correct trials (non-significant ones are not accounted for)
    MminDe(c) = nanmean(minD{c}(~k));                                       % and for error trials
    
    FminDc(c) = 100*sum(~isnan(minD{c}(k))) / sum(k);                       % Count ratio of significant hyperpolarizations over correct trials
    FminDe(c) = 100*sum(~isnan(minD{c}(~k))) / sum(~k);                     % and over error trials
end

%% HYPERPOL AMPLITUDE IN CORRECT VS ERROR
figure;
subplot(221); hold on; 
for c = 1:Nr
    plot(1:2,[MminDc(c) MminDe(c)],'o-','Color',[0.8 0.8 0.8]);
end
[~,p_hyper] = ttest(MminDc,MminDe)
errorbar([0.9 2.1],[nanmean(MminDc) nanmean(MminDe)],[SEM(MminDc) SEM(MminDe)],'ok-')
plot_significance(p_hyper,0.9,2.1,nanmean(MminDc),nanmean(MminDe))
xlim([0.8 2.2]);
ylabel('Hyperpol ampl');

%% HYPERPOL FREQUENCY IN CORRECT VS ERROR
subplot(223); hold on; 
for c = 1:Nr
    plot(1:2,[FminDc(c) FminDe(c)],'o-','Color',[0.8 0.8 0.8]);
end
[~,p_freq] = ttest(FminDc,FminDe)
errorbar([0.9 2.1],[nanmean(FminDc) nanmean(FminDe)],[SEM(FminDc) SEM(FminDe)],'ok-')
plot_significance(p_freq,0.9,2.1,nanmean(FminDc),nanmean(FminDe))
xlim([0.8 2.2]);
ylabel('Hyperpol freq');

%% PLOT HYPERPOLARIZATION AMPLITUDE VS PERFORMANCE
subplot(222); hold on; 
scatter(Perf, MminD,'ok')

% LEAST SQUARES ESTIMATE FOR DISTRIBUTION
k = ~isnan(MminD);                                                          % REMOVE CELLS WITH NO HYPERPOL
coefficients = polyfit(Perf(k), MminD(k), 1);                               % Equivalent to: [ones(size(D)), D]\C
Cfitted = polyval(coefficients, [min(Perf(k)) max(Perf(k))]);
pfit = coefTest(fitlm(Perf(k), MminD(k),'linear'));
plot([min(Perf(k)) max(Perf(k))] , Cfitted,'-k')
title(num2str(pfit));

% mdl = fitlm(Perf,MminD)

%% PLOT HYPERPOLARIZATION FREQUENCY VS PERFORMANCE
subplot(224); hold on; 
scatter(Perf, FminD,'ok')

% LEAST SQUARES ESTIMATE FOR DISTRIBUTION
coefficients = polyfit(Perf, FminD, 1);                                     % Equivalent to: [ones(size(D)), D]\C
Cfitted = polyval(coefficients, [min(Perf) max(Perf)]);
pfit = coefTest(fitlm(Perf, FminD,'linear'));
plot([min(Perf) max(Perf)] , Cfitted,'-k')
title(num2str(pfit));

% mdl = fitlm(Perf,FminD)

