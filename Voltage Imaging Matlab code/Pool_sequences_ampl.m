function Pool_sequences_ampl(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
elseif strcmp(INclass,'SOM')
    Days_SOM;
end

cols = cell_colors;
freqtit = {'Delta (0.5-3Hz)';'Theta (4-10Hz)';'Beta (15-30Hz)';'Gamma (40-90Hz)'};

MA = cell(4,3);
% MF0 = cell(1,3);
% MFall = cell(1,3);
% MFn = cell(1,3);
Fp = cell(1,2);

Fs = 1000;
time = 0 : 1/Fs : (11*Fs-1)/Fs;

%% LOAD POOLED DATA AND STACK IN SINGLE COLUMN
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'AMP','MC','TR','timepoints');

AMP = vertcat(AMP{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

ls = length(AMP);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
AMPall = {};
MCall = [];
TRall = {};

for i = 1:ls
    for r = 1:size(AMP{i},1)
        AMPall{count,1} = permute(AMP{i}(r,:,:,:),[3,2,4,1]);                          % Keep only that cell's DFF [trials x time]
        MCall(count,:) = MC{i}(r,:);    
        TRall{count,1} = TR{i};
         count = count + 1;
    end
end
AMP = AMPall;
MC = MCall;
TR = TRall;
clear AMPall MCall TRall

ls = count-1;

%% ZSCORE 
ontime = (time >= timepoints(1) & time < timepoints(3));
for i = 1:ls   
    AMP{i} = (AMP{i} - mean(AMP{i}(:,ontime,:),2)) ./ std(AMP{i}(:,ontime,:),[],2);     % ZSCORE OVER MODULATION BINS ONLY!!!
    AMP{i}(isnan(AMP{i}) | isinf(AMP{i})) = 0;
end
          
%% KEEP MEAN AMPLITUDE OVER EACH ODOR AND COLLECTIVELY FOR EACH CELL GROUP
for tt1 = 1:4                                                               % For each cell group
    k = (MC(:,1) == tt1);                                                   % Keep corresponding cells
    if tt1 == 3
        k = (MC(:,1) == 0);
    elseif tt1 == 4
        k = isnan(MC(:,1));
    end
    
    amp = AMP(k);                                                           % Keep their DFF
    tr = TR(k);
    
    for i = 1:length(amp)                                                   % For each cell
        for tt2 = 1:2                                                       % For each first odor
            ma = amp{i}(tr{i} == tt2,:,:);                                    % Keep those trials and all freqs
            ma = mean(ma,1);                                                % Average over trials  [1 x time]
            MA{tt1,tt2} = [MA{tt1,tt2}; ma];                                % Stack in corresponding cell array
        end
        ma = mean(amp{i},1);                                                % Same for ALL TRIALS
        MA{tt1,3} = [MA{tt1,3}; ma];
    end
end

%% SORT FIELDS
for tt1 = 1:4                                                               % For each cell group
    k = (MC(:,1) == tt1);                                                   % Keep corresponding cells
    if tt1 == 3
        k = (MC(:,1) == 0);
    elseif tt1 == 4
        k = isnan(MC(:,1));
    end
    [Fp{tt1},order] = sort(MC(k,2));                                        % Sort their fields
        
    for tt2 = 1:3
        MA{tt1,tt2} = MA{tt1,tt2}(order,:,:);
    end
end

%% PLOT AVERAGE AMPLITUDES OF EACH GROUP OVER EACH ODOR AND IN ALL TRIALS
ylabels = {'OdorA';'OdorB';'Non-spec';'Non-field'};
tits = {'OdorA','OdorB','All trials'};

for f = 1:4
    figure('Name',['Average Amplitude - ',freqtit{f}]);
    for tt1 = 1:4                                                  % For each trial type
        for tt2 = 1:3                                              % For each trial type (again)
            subplot(4,3,(tt1-1)*3 + tt2); hold on;
            plot_seq_rates(MA{tt1,tt2}(:,:,f),[],time,timepoints,1); % Plot pooled sequences
            plot(Fp{tt1},1:length(Fp{tt1}),'.k')
            
            xlim([0, timepoints(5)]);
            ylabel(ylabels{tt1});
            if tt1 == 1, title(tits{tt2}); end
        end
    end
end

%% COMPARE PREFERRED AND NON-PREFFERED TRIALS OF ODOR-SPECIFIC MCELLS
MApref = [MA{1,1}; MA{2,2}];
MAnpref = [MA{1,2}; MA{2,1}];

% SORT BY FIELDS
Fmc = [Fp{1};Fp{2}];
[Fmc,order] = sort(Fmc);
MApref = MApref(order,:,:);
MAnpref = MAnpref(order,:,:);
%------
    
for f = 1:4
    fig(f) = figure('Name',freqtit{f});
    
    subplot(4,4,[1 5 9]);
    plot_seq_rates(MApref(:,:,f),[],time,timepoints,1);
    plot(Fmc,1:length(Fmc),'.k')
    xlim([time(1), timepoints(5)]);
    title('Preferred odor');
    ylabel('All odor-cells');
    
    subplot(4,4,[2 6 10]);
    plot_seq_rates(MAnpref(:,:,f),[],time,timepoints,1);
    xlim([time(1), timepoints(5)]);
    title('Non-Preferred odor');
    
    subplot(4,4,[13 14]); hold on;
    plot_mrate_traces(nanmean(MApref(:,:,f),1),SEM(MApref(:,:,f)),time,timepoints,'k');
    plot_mrate_traces(nanmean(MAnpref(:,:,f),1),SEM(MAnpref(:,:,f)),time,timepoints,'r');
    
    % SIGNIFICANCE DURING ODOR 
    ontime = find(time >= timepoints(1) & time < timepoints(2));
    p_pnp = nan(1,length(ontime));
    for i = 1:length(ontime)
        p_pnp(i) = ranksum(MApref(:,ontime(i),f),MAnpref(:,ontime(i),f));
    end
    [~,p_pnp] = fdr(p_pnp);
    
    tontime = time(ontime);
    plot(tontime(p_pnp<0.05),1.1*max(mean(MApref(:,ontime,f),1))*ones(sum(p_pnp<0.05),1),'*k')
end

%% COMPARE POOLED MCELL AND NON-MCELL AMPLITUDES
MAmc = cell2mat(MA(1:3,3));                                                 % Pool ALL odor-cells at ALL trials
MAnmc = MA{4,3};                                                            % And keep nonMcells at ALL trials

% SORT ODORCELLS BY FIELDS
Fmc = [Fp{1};Fp{2};Fp{3}];
[Fmc,order] = sort(Fmc);
MAmc = MAmc(order,:,:);
%------

for f = 1:4
    figure(fig(f));
    
    subplot(4,4,[3 7 11]);
    plot_seq_rates(MAmc(:,:,f),[],time,timepoints,1);
    plot(Fmc,1:length(Fmc),'.k')
    xlim([time(1), timepoints(5)]);
    title('Mcells');
    
    subplot(4,4,[4 8 12]);
    plot_seq_rates(MAnmc(:,:,f),[],time,timepoints,1);
    plot(Fp{4},1:length(Fp{4}),'.k')
    xlim([time(1), timepoints(5)]);
    title('Non Mcells');
    
    subplot(4,4,[15 16]); hold on;
    plot_mrate_traces(nanmean(MAmc(:,:,f),1),SEM(MAmc(:,:,f)),time,timepoints,'k');
    plot_mrate_traces(nanmean(MAnmc(:,:,f),1),SEM(MAnmc(:,:,f)),time,timepoints,'r');
    
    % SIGNIFICANCE DURING ODOR 
    ontime = find(time >= timepoints(1) & time < timepoints(2));
    p_amp = nan(1,length(ontime));
    for i = 1:length(ontime)
        p_amp(i) = ranksum(MAmc(:,ontime(i),f),MAnmc(:,ontime(i),f));
    end
    [~,p_amp] = fdr(p_amp);
    
    tontime = time(ontime);
    plot(tontime(p_amp<0.05),1.1*max(mean(MAmc(:,ontime,f),1))*ones(sum(p_amp<0.05),1),'*k')
end


%% PLOT ALL CELLS IN ALL TRIALS (REVISIONS)
MA = cell2mat(MA(:,3));                                                 % Pool ALL cells at ALL trials

figure;
for f = 1:4
   subplot(4,1,f); hold on;
   plot_mrate_traces(nanmean(MA(:,:,f),1),SEM(MA(:,:,f)),time,timepoints,'k');
end


