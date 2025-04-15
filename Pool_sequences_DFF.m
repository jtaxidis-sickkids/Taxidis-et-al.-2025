function [MD,FMD,Fp] = Pool_sequences_DFF(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
elseif strcmp(INclass,'SOM')
    Days_SOM;
end

MD = cell(4,3);
FMD = cell(4,3);
Fp = cell(4,1);

freq = [0.5 3;        % Delta (based on power spectra early peak)
        4   10;       % Theta frequency range (Bittner et al 2015 but also the power spectra peaks).
        15  30;      % Beta (based on literature and mean ISI distributions)
        40  90];     % Gamma

Fs = 1000;
for f = 1:4
    [filtb{f},filta{f}] = butter(3,freq(f,:)/(Fs/2),'bandpass');                        % construct the bandpass filter
end
freqtit = {'Delta (0.5-3Hz)';'Theta (4-10Hz)';'Beta (15-30Hz)';'Gamma (40-90Hz)'};

%% LOAD POOLED DATA AND STACK
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'DFF','MC','TR','timepoints','time');
    
DFF = vertcat(DFF{:});
MC = vertcat(MC{:});    
TR = vertcat(TR{:});

ls = length(MC);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
DFFall = {};
MCall = [];
TRall = {};

for i = 1:ls
    for r = 1:size(DFF{i},1)
        DFFall{count,1} = squeeze(DFF{i}(r,:,:))';                          % Keep only that cell's DFF [trials x time]
        MCall(count,:) = MC{i}(r,:);    
        TRall{count,1} = TR{i};
        count = count + 1;
    end
end
DFF = DFFall;
MC = MCall;
TR = TRall;
clear DFFall MCall TRall

ls = count-1;

%% SMOOTH RAW DFF
for i = 1:ls
    DFF{i} = sgolayfilt(DFF{i}, 1, 21, [] , 2);                            % Apply linear filter along time bins
end

%% FILTER DFF
FDFF = cell(size(DFF));
for i = 1:ls
    FDFF{i} = nan([size(DFF{i}),4]);
    for tr = 1:size(DFF{i},1)
        for f = 1:4
            FDFF{i}(tr,:,f) = filtfilt(filtb{f},filta{f},DFF{i}(tr,:));     % Apply linear filter along time bins
        end
    end
end

%% KEEP MEAN DFF AND FILTERED DFF OF ALL CELLS AND TRIALS
DFFall = cellfun(@(x) mean(x,1), DFF,'UniformOutput',false);                % Mean per cell across all trials
DFFall = cell2mat(DFFall);

FDFFall = cellfun(@(x) mean(x,1), FDFF,'UniformOutput',false);              % Same for each freq. separately
FDFFall = cell2mat(FDFFall);

%% KEEP MEAN DFF OVER EACH ODOR AND COLLECTIVELY FOR EACH CELL GROUP
for tt1 = 1:4                                                               % For each cell group
    k = (MC(:,1) == tt1);                                                   % Keep corresponding cells
    if tt1 == 3
        k = (MC(:,1) == 0);
    elseif tt1 == 4
        k = isnan(MC(:,1));
    end
    
    dff = DFF(k);                                                           % Keep their DFF
    tr = TR(k);
    
    for i = 1:length(dff)                                                   % For each cell
        for tt2 = 1:2                                                       % For each first odor
            md = dff{i}(tr{i} == tt2,:);                                    % Keep those trials
            md = mean(md,1);                                                % Average over trials  [1 x time]
            MD{tt1,tt2} = [MD{tt1,tt2}; md];                                % Stack in corresponding cell array
        end
        md = mean(dff{i},1);                                                % Same for ALL TRIALS
        MD{tt1,3} = [MD{tt1,3}; md];
    end
    
    
    % REPEAT FOR FILTERED DFF
    dff = FDFF(k);                                                           % Keep their DFF
    for i = 1:length(dff)                                                   % For each cell
        for tt2 = 1:2                                                       % For each first odor
            md = dff{i}(tr{i} == tt2,:,:);                                    % Keep those trials
            md = mean(md,1);                                                % Average over trials  [1 x time]
            FMD{tt1,tt2} = [FMD{tt1,tt2}; md];                                % Stack in corresponding cell array
        end
        md = mean(dff{i},1);                                                % Same for ALL TRIALS
        FMD{tt1,3} = [FMD{tt1,3}; md];
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
        
    for tt2 = 1:3                                                           % For each first odor
        MD{tt1,tt2} = MD{tt1,tt2}(order,:);
        
        FMD{tt1,tt2} = FMD{tt1,tt2}(order,:,:);
    end
end

%% PLOT AVERAGE DFF OF EACH GROUP OVER EACH ODOR AND IN ALL TRIALS
figure('Name','Average DFFs');
ylabels = {'OdorA';'OdorB';'Non-spec';'Non-field'};
tits = {'OdorA','OdorB','All trials'};
for tt1 = 1:4                                                  % For each trial type
    for tt2 = 1:3                                              % For each trial type (again)
        subplot(4,3,(tt1-1)*3 + tt2); hold on;
        plot_seq_rates(MD{tt1,tt2},[],time,timepoints,0.2); % Plot pooled sequences
        plot(Fp{tt1},1:length(Fp{tt1}),'.k')
        
        xlim([0, timepoints(5)]);
        ylabel(ylabels{tt1});
        if tt1 == 1, title(tits{tt2}); end
    end
end

%% PLOT MEAN DFF AND FILTERED DFF OF ALL CELLS OVER ALL TRIALS
figure;
subplot(511);
plot_mrate_traces(nanmean(DFFall,1),SEM(DFFall),time,timepoints,'k');       % Plot mean across cells
ylabel('Mean DFF ALL cells');

for f = 1:4
   subplot(5,1,f+1);
    plot_mrate_traces(nanmean(FDFFall(:,:,f),1),SEM(FDFFall(:,:,f)),time,timepoints,'k'); % Plot mean across cells for each freq.
    ylabel(['Mean ',freqtit{f},' ALL cells']);   
end

%% COMPARE PREFERRED AND NON-PREFFERED TRIALS OF ODOR-SPECIFIC MCELLS
figure;

MDpref = [MD{1,1}; MD{2,2}];
MDnpref = [MD{1,2}; MD{2,1}];

% SORT BY FIELDS
Fmc = [Fp{1};Fp{2}];
[Fmc,order] = sort(Fmc);
MDpref = MDpref(order,:);
MDnpref = MDnpref(order,:);
%------

subplot(4,4,[1 5 9]);
plot_seq_rates(MDpref,[],time,timepoints,0.2);     
plot(Fmc,1:length(Fmc),'.k')
xlim([time(1), timepoints(5)]);
title('Preferred odor');
ylabel('All odor-cells');

subplot(4,4,[2 6 10]);
plot_seq_rates(MDnpref,[],time,timepoints,0.2);     
xlim([time(1), timepoints(5)]);
title('Non-Preferred odor');

subplot(4,4,[13 14]); hold on;
plot_mrate_traces(nanmean(MDpref,1),SEM(MDpref),time,timepoints,'k');
plot_mrate_traces(nanmean(MDnpref,1),SEM(MDnpref),time,timepoints,'r');

% SIGNIFICANCE DURING ODOR + 1 sec DELAY
ontime = find(time >= timepoints(1) & time < timepoints(2) + 1);
p_pnp = nan(1,length(ontime));
for i = 1:length(ontime)
    p_pnp(i) = ranksum(MDpref(:,ontime(i)),MDnpref(:,ontime(i)));
end
[~,p_pnp] = fdr(p_pnp);

tontime = time(ontime);
plot(tontime(p_pnp<0.05),1.1*max(mean(MDpref,1))*ones(sum(p_pnp<0.05),1),'*k')

%% COMPARE POOLED MCELL AND NON-MCELL DFF
MDmc = cell2mat(MD(1:3,3));                                                 % Pool ALL odor-cells at ALL trials
MDnmc = MD{4,3};                                                            % And keep nonMcells at ALL trials

% SORT ODORCELLS BY FIELDS
Fmc = [Fp{1};Fp{2};Fp{3}];
[Fmc,order] = sort(Fmc);
MDmc = MDmc(order,:);
%------

subplot(4,4,[3 7 11]);
plot_seq_rates(MDmc,[],time,timepoints,0.2);     
plot(Fmc,1:length(Fmc),'.k')
xlim([time(1), timepoints(5)]);
title('Mcells');

subplot(4,4,[4 8 12]);
plot_seq_rates(MDnmc,[],time,timepoints,0.2); 
plot(Fp{4},1:length(Fp{4}),'.k')
xlim([time(1), timepoints(5)]);
title('Non Mcells');

subplot(4,4,[15 16]); hold on;
plot_mrate_traces(nanmean(MDmc,1),SEM(MDmc),time,timepoints,'k');
plot_mrate_traces(nanmean(MDnmc,1),SEM(MDnmc),time,timepoints,'r');

% SIGNIFICANCE DURING ODOR + 1 sec DELAY
ontime = find(time >= timepoints(1) & time < timepoints(2) + 1);
p_dff = nan(1,length(ontime));
for i = 1:length(ontime)
    p_dff(i) = ranksum(MDmc(:,ontime(i)),MDnmc(:,ontime(i)));
end
[~,p_dff] = fdr(p_dff);

tontime = time(ontime);
plot(tontime(p_dff<0.05),1.1*max(mean(MDmc,1))*ones(sum(p_dff<0.05),1),'*k')


%% COMPARE DELTA AND THETA FILTERED DFF OF POOLED MCELL AND NON-MCELL
FMDmc = cell2mat(FMD(1:3,3));                                                 % Pool ALL odor-cells at ALL trials
FMDnmc = FMD{4,3};                                                            % And keep nonMcells at ALL trials
FMDmc = FMDmc(order,:,:);

figure;
for f = 1:4
    Dmc = FMDmc(:,:,f);
    Dnmc = FMDnmc(:,:,f);
    subplot(4,8,[1 9 17]+2*f-2);
    plot_seq_rates(Dmc,[],time,timepoints,0.2);
    plot(Fmc,1:length(Fmc),'.k')
    xlim([time(1), timepoints(5)]);
    title('Odor-specific Mcells');
    ylabel(freqtit{f});
    
    subplot(4,8,[1 9 17]+2*f-1);
    plot_seq_rates(Dnmc,[],time,timepoints,0.2);
    plot(Fp{4},1:length(Fp{4}),'.k')
    xlim([time(1), timepoints(5)]);
    title('Non Mcells');
    
    subplot(4,8,[25 26]+2*f-2); hold on;
    plot_mrate_traces(nanmean(Dmc,1),SEM(Dmc),time,timepoints,'k');
    plot_mrate_traces(nanmean(Dnmc,1),SEM(Dnmc),time,timepoints,'r');
    
    % SIGNIFICANCE DURING ODOR + 1 sec DELAY
    ontime = find(time >= timepoints(1) & time < timepoints(2) + 1);
    p_dff = nan(1,length(ontime));
    for i = 1:length(ontime)
        p_dff(i) = ranksum(Dmc(:,ontime(i)),Dnmc(:,ontime(i)));
    end
    [~,p_dff] = fdr(p_dff);
    
    tontime = time(ontime);
    plot(tontime(p_dff<0.05),1.1*max(mean(Dmc,1))*ones(sum(p_dff<0.05),1),'*k')
end


