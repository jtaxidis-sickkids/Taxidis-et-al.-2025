%% PLOT FOV
Days_Step3_pooled
for i = [2 7 30]
    day = s3days{i};                                                         % Plot part of FOV of JT88 on 2017_12_27
    
    lims = [30 482];
    
    figure;
    subplot(121);hold on;
    plot_FOV(day,lims);                                                         % Plot the 2-channel FOV
    line(10 + [0 102.4],[460 460],'color','w','linewidth',4)                     % Draw line scale of 100µm (102.4 pixels)
    
    subplot(122);hold on;
    plot_segmented_FOV(day,lims);
end

%% ORGANIZE AND SAVE THE DECONVOLVED SIGNAL IN TRIALS STRUCTURE
Days_Step3_pooled;
for d = 1:length(s3days)
    [animal,main_folder,day] = get_animal_name(s3days{d});
    s3days{d} = fullfile(main_folder,'Processed_files',day);
end
Organize_SpikeSignal(s3days);

%% DETECT MCELLS BASED ON DECONVOLVED SIGNAL
Days_Step3_pooled;
s3days([1 7 13 32 40 48 53 54]) = [];   % Remove days requiring multiple Modulation detections

for d = 1:length(s3days)
    Get_Modulation_SpikeSignal(s3days{d})
end

%% PLOT POOLED SEQUENCES
Days_Step3_pooled;
s3days([1 7 13 32 40 48 53 54]) = [];   % Remove days requiring multiple Modulation detections

[MP,MNP,FF] = Pool_sequences_SpikeSignal(s3days,'firstodor');

%%  CELLS SORTED BY MAXIMUM MEAN SIGNAL OVER THE TWO ODORS
Days_Step3_pooled;
s3days([1 7 13 32 40 48 53 54]) = [];   % Remove days requiring multiple Modulation detections

load(fullfile(s3days{1},'Modulation_delay_5_firstodor_ASAP.mat')); % Load info for all sessions of that delay split by that criterion
load(fullfile(s3days{5},'Trials_delay_5.mat'),'time');

tend = find(time > timepoints(6),1,'first');
time = time(1:tend);
odor1 = (time > timepoints(1) & time < timepoints(2));                      % Get timepoints of the first odor

MS = cell(1,2);
SS = cell(1,2);
zMS = cell(1,2);
MC = cell(1,2);
SIn = cell(1,2);
Nr = zeros(2,1);
zMODOR = cell(1,2);

for ct = 1:2
    MS{ct} = [];
    SS{ct} = [];
    MC{ct} = [];
    SIn{ct} = [];
    
    for d = 1:length(s3days)                                                % For each day
        load(fullfile(s3days{d},'Trials_delay_5.mat'),'Npy');      % Load necessary data
        load(fullfile(s3days{d},'Trials_delay_5_SpikeSignal.mat'),'S');
        load(fullfile(s3days{d},'Modulation_delay_5_firstodor_ASAP.mat'),'Mcells','SI');
        
        S = R_ctype(S,ct,Npy);                                              % Keep the signal of corresponding cell type
        ms = mean(S(:,1:tend,:),3);                                         % Get mean signal over all trials
        ss = std(S(:,1:tend,:),[],3) / sqrt(size(S,3));                     % Get SE signal over all trials

        if size(Mcells,1) == 1, Mcells{2,1} = zeros(0,3);Mcells{2,2} = zeros(0,3);end
        mc = Mcells(ct,:);                                                  % Keep the Mcells of that cell type
        
        f = nan(size(S,1),2);                                               % Create fields matrix for all cells
        si = nan(size(S,1),1);                                              % and SI matrix
        for tt = 1:2                                                        % For each first odor
            if isempty(mc{tt}), mc{tt} = zeros(0,3);end
            
            f(mc{tt}(:,1),1) = tt;                                          % Tag that odor only in the indexes of mcells
            f(mc{tt}(:,1),2) = mc{tt}(:,2);                                 % and tag the second column with their fields
            
            si(mc{tt}(:,1)) = SI{ct,tt};                                    % Store the SI of Mcells
        end
        
        MS{ct} = [MS{ct}; ms];                                              % Store mean signal traces
        SS{ct} = [SS{ct}; ss];                                              % Store mean signal traces
        MC{ct} = [MC{ct}; f];                                               % Store fields
        SIn{ct} = [SIn{ct}; si];                                            % Store SI
    end
    
    zMS{ct} = zscore(MS{ct},[],2);                                           % ZSCORE THE MEAN SIGNALS
    
    % SORT BY (NON ZSCORED) MEAN SIGNALS DURING ODOR
    modor = mean(MS{ct}(:,odor1),2);                                        % Mean signal over odor bins
    
    [~,order] = sort(modor,'descend');                                      % Sort the mean signals over odors (descending order)
    MS{ct} = MS{ct}(order,:);                                               % Reorder accordingly
    SS{ct} = SS{ct}(order,:);                                               % Reorder accordingly
    zMS{ct} = zMS{ct}(order,:);                                             % Reorder accordingly
    MC{ct} = MC{ct}(order,:);
    SIn{ct} = SIn{ct}(order);
    % -----------------
    
    zMODOR{ct} = mean(zMS{ct}(:,odor1),2);                                  % Mean zscored signal over the odor
    
    Nr(ct) = size(MS{ct},1);
end

%% COUNT TIME CELLS WITH NEGATIVE ODOR RESPONSES
for ct = 1:2
    z(ct) = 100 * sum(zMODOR{ct} < 0) / length(zMODOR{ct});                              % Compute % of all cells with negative odor response
    
    f = MC{ct}(:,2);                                                        % All fields from all Mcells
    dmc = find(f > timepoints(2));                                          % Find indexes of delay cells
    zz = zMODOR{ct}(dmc);                                                   % Keep their relative odor responses  
    z_de(ct) = 100 * sum(zz < 0) / length(zz);                           % Compute % of time cells with negative odor response
end
z
z_de

%% PLOT
figure;
C = cell_colors;
Ctypes = {'PY';'IN'};

for ct = 1:2                                                                % For each cell type
    % PLOT MEAN ZSCORED SIGNALS
    subplot(5,8,(ct-1)*4 + [9:11 17:19 25:27]); hold on;
    plot_seq_rates(zMS{ct},[],[],ontime,time,timepoints,2);
    for tt = 1:2
        k = find(MC{ct}(:,1) == tt);
        plot(MC{ct}(k,2),k,'ok','MarkerFaceColor',C(tt,:),'Markersize',3);
    end
    xlim([0,timepoints(6)]);
    ylabel(Ctypes{ct});
    
    % PLOT MEAN ZSCORED SIGNALS ACROSS CELLS
    subplot(5,8,(ct-1)*4 + [1:3]); hold on;
    fill_plot(time,mean(zMS{ct},1),std(zMS{ct},[],1)/sqrt(Nr(ct)),'k');
    axis tight;
    xlim([0,timepoints(6)]);
    
    % PLOT MEAN ZSCORED SIGNALS ACROSS CELLS WITH POSITIVE/NEGATIVE ODOR RESPONSE
    subplot(5,8,(ct-1)*4 + [1:3]); hold on;
    k = zMODOR{ct}>0;
    fill_plot(time,mean(zMS{ct}(k,:),1),std(zMS{ct}(k,:),[],1)/sqrt(sum(k)),'b');
    k = zMODOR{ct}<=0;
    fill_plot(time,mean(zMS{ct}(k,:),1),std(zMS{ct}(k,:),[],1)/sqrt(sum(k)),'r');
    axis tight;
    xlim([0,timepoints(6)]);
    
    % PLOT MEAN SIGNALS ACROSS CELLS ZOOMED AROUND ODOR
    subplot(5,8,(ct-1)*4 + 4); hold on;
    k1 = find(time >= 0.5,1,'first');
    k2 = find(time <= 2.5,1,'last');
    fill_plot(time(k1:k2),mean(zMS{ct}(:,k1:k2),1),std(zMS{ct}(:,k1:k2),[],1)/sqrt(Nr(ct)),'k');
    k = zMODOR{ct}>0;
    fill_plot(time(k1:k2),mean(zMS{ct}(k,k1:k2),1),std(zMS{ct}(k,k1:k2),[],1)/sqrt(sum(k)),'b');
    k = zMODOR{ct}<=0;
    fill_plot(time(k1:k2),mean(zMS{ct}(k,k1:k2),1),std(zMS{ct}(k,k1:k2),[],1)/sqrt(sum(k)),'r');
    axis tight;
    xlim([time(k1), time(k2)]);
    set(gca,'XTick',0.6:0.2:2.4);
    
    % PLOT MEAN ZSCORED SIGNAL OVER THE ODOR
    subplot(5,8,(ct-1)*4 + [12 20 28]); hold on; view(90, 90);
    plot(1:Nr(ct),zMODOR{ct},'.','Markersize',4);
    line([1 Nr(ct)],[0 0],'LineStyle','--','Color','k');
    axis tight;
    
    edges = 1:200:Nr(ct);                                                   % Make bins
    de = edges(2)-edges(1);                                                 % bin size
    [M,S] = bin_x_axis(1:Nr(ct),zMODOR{ct},edges);                              % Bin the MEAN ZSCORED SIGNAL
    fill_plot(edges(1:end-1)+de/2, M, S, 'b');                              % Plot
    
    % PLOT SI
    subplot(5,8,(ct-1)*4 + [33:35]); hold on;
    f = MC{ct}(~isnan(SIn{ct}),2);                                          % Keep all Mcell fields
    si = SIn{ct}(~isnan(SIn{ct}));                                          % and their SI
    plot(f + randn(size(f))*0.06 , si + randn(size(si))*0.01,'.','Markersize',4); % Plot SI as function of field (with jitter)
    
    ton = time(ontime);                                                     % Keep ontime only
    edges = ton(1:3:end);                                                   % Make bins
    de = edges(2)-edges(1);                                                 % bin size
    [M,S] = bin_x_axis(f,si,edges);                                         % Bin the SI distribution
    fill_plot(edges(1:end-1)+de/2, M, S, 'b');                              % Plot
    ylim([0,1.2]); xlim([0.9,7.1]);
    
    % COMPARE SI BETWEEN ODOR ONSET AND REST
    thresh  = 0.2;

    subplot(5,8,(ct-1)*4 + 36); hold on;
    si1 = si(f > timepoints(1) & f <= timepoints(1) + thresh);
    si2 = si(f > timepoints(1) + thresh);
    p = ranksum(si1,si2,'tail','left');                                  % left-tailed two-sample ttest
    
    dp = table(si1',si2','VariableNames',{'Onset','Post'});
    violinplot(dp);
    hold on;
    plot_significance(p,1,2,si1,si2);
    title(num2str(p));
end

%% PLOT EXAMPLE SPIKE SIGNALS (TOP AND BOTTOM  ODOR RESPONSES)
figure;

fMS = flipud(MS{1});                                                        % Flip the mean signal matrix to place the highest responses on top
fSS = flipud(SS{1});                                                        % Flip the mean signal matrix to place the highest responses on top
cell_list = [2.1 2.4 2.8 2.9 3.2 3.3 3.6 4.5,...                            % List of example odor-negative response cells
    10.7 10.9 11.3 11.4 11.6 12.1 12.2 12.3] * 1e3;                         % List of example odor-positive response cells

shade_timepoints([0 30],timepoints);                                        % Make shaded rectangles for timepoints
for c = 1:16                                                                % For each cell
    fms = smooth(fMS(cell_list(c),:),5,'moving');                           % SMOOTH THE MEAN SIGNAL
    fss = smooth(fSS(cell_list(c),:),5,'moving');                           % SMOOTH THE SE OF SIGNAL
    fill_plot(time,fms'*100 + c,fss'*100,'k');                                     % Plot shaded rate
end
axis tight
xlim([0 9]);