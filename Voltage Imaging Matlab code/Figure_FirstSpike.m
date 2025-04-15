function Figure_FirstSpike(INclass,t_firstspike)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
elseif strcmp(INclass,'SOM')
    Days_SOM;
end

%% LOAD POOLED DATA AND STACK IN SINGLE COLUMNS
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'SP','DFF','MC','TR','time');

SP = vertcat(SP{:});
DFF = vertcat(DFF{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

ls = length(SP); % Number of sessions

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
SPall = {};
DFFall = {};
MCall = [];
TRall = {};

for i = 1:ls
    for r = 1:size(SP{i},1)
        DFFall{count,1} = squeeze(DFF{i}(r,:,:))';                          % Keep only that cell's DFF [trials x time]
        SPall{count,1} = SP{i}(r,:);
        count = count + 1;
    end
end
DFF = DFFall;
SP = SPall;
clear SPall MCall TRall DFFall

ls = count-1;

%% COUNT FIRST-SPIKES PER TRIAL
FSP = cell(ls,1);
for i = 1:ls                                                                % For each cell
    FSP{i} = nan(length(SP{i}),1);
    for tr = 1:length(SP{i})                                                % For each trial
        sp = SP{i}{tr} ;                                                    % Keep all spikes
        firstspike = (sp >= t_firstspike(1) & sp <= t_firstspike(2));       % Find spikes within the window
        FSP{i}(tr) = sum(firstspike);                                       % Count them
    end
end

%% FIND REBOUND ONSET TIME PER TRIAL
RB = cell(ls,1);
for i = 1:ls                                                                % For each cell
    RB{i} = nan(length(SP{i}),1);
    for tr = 1:length(SP{i})                                                % For each trial
        sp = SP{i}{tr} ;                                                    % Keep all spikes
        k = find(sp > t_firstspike(2) & sp < 2);
        if ~isempty(k)
            RB{i}(tr) = (sp(k(1)));                   %
        end
    end
end

%% GROUP DFF AND SPIKES INTO 2 CATEGORIES
DFFon = cell(ls,1);
DFFoff = cell(ls,1);

Ron = cell(ls,1);
Roff = cell(ls,1);

for i = 1:ls
    DFFon{i} = DFF{i}(FSP{i} > 0,:);                                        % DFFs in trials with first spike
    DFFoff{i} = DFF{i}(FSP{i} == 0,:);                                      % DFFs in trials with no first spike

    SPon = SP{i}(FSP{i} > 0);
    SPoff = SP{i}(FSP{i} == 0);
    
    for tr = 1:sum(FSP{i} > 0)
        Ron{i}(tr,:) = histcounts(SPon{tr},0:0.005:10) / 0.005;             % Firing rates (5msec bin) in trials with first spike
    end
    for tr = 1:sum(FSP{i} == 0)
        Roff{i}(tr,:) = histcounts(SPoff{tr},0:0.005:10) / 0.005;
    end
end
Ron = cellfun(@(x) nanmean(x,1), Ron,'UniformOutput',false);                   % MEan rates over trials
Roff = cellfun(@(x) nanmean(x,1), Roff,'UniformOutput',false);

DFFon = cellfun(@(x) nanmean(x,1), DFFon,'UniformOutput',false);                   % MEan rates over trials
DFFoff = cellfun(@(x) nanmean(x,1), DFFoff,'UniformOutput',false);

kon = cellfun(@isempty,Ron);
koff = cellfun(@isempty,Roff);

Ron(kon) = [];
Roff(koff) = [];
DFFon(kon) = [];
DFFoff(koff) = [];

Ron = cell2mat(Ron);
Roff = cell2mat(Roff);
DFFon = cell2mat(DFFon);
DFFoff = cell2mat(DFFoff);

%% PLOT MEAN DFF WITH VS WITHOUT FIRST SPIKE
figure;

subplot(311); hold on;
histogram('BinEdges',0:0.005:10,'BinCounts',mean(Ron,1));
histogram('BinEdges',0:0.005:10,'BinCounts',mean(Roff,1));
xlim([0.8 2.2]);

subplot(312); hold on;
fill_plot(time,mean(DFFon),SEM(DFFon),'k');
fill_plot(time,mean(DFFoff),SEM(DFFoff),'r');
xlim([0.8, 2.2]);

%% PLOT RELATIONSHIP BETWEEN FIRST SPIKE AND TIME OF ONSET SPIKING
fsp = cell2mat(FSP);
rb = cell2mat(RB);

subplot(313); hold on  
scatter(fsp+randn(size(fsp))*0.1,rb,'.k');
[M,S,D] = bin_x_axis(fsp,rb,0:1:max(fsp)+1);
fill_plot(0:1:max(fsp),M,S,'r')
xlabel('amount of first spike');
ylabel('time of rebound spiking')

