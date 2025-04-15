function [MR, order] = Plot_Rates_pooled(sessions,optnum,optoON,order)

%% SET UP PARAMETERS
ls = length(sessions);

tp = [1 2 7 8 9 11 11];
Ctypes = {'PY';'IN'};
C = cell_colors;

if isempty(order)
    order = cell(2,1);
end

%%  CELLS SORTED BY MAXIMUM MEAN SIGNAL OVER THE TWO ODORS
dirname = ['../ProcessedData/',sessions{1}];
load(fullfile(dirname,'CA1data.mat'),'bins');                          
trialbins = (bins < 11);
bins = bins(trialbins);
onbins = (bins > 1 & bins < 7);
odor1 = (bins > 1 & bins < 2);                                              % Get timepoints of the first odor
baseline = (bins < 1);

MR = cell(1,2);
SR = cell(1,2);
MC = cell(1,2);
SIn = cell(1,2);
Nr = zeros(2,1);
zMODOR = cell(1,2);

for ct = 1:2
    MR{ct} = [];
    SR{ct} = [];
    MC{ct} = [];
    SIn{ct} = [];
    
    for s = 1:ls                                                            % For each session
        dirname = ['../ProcessedData/',sessions{s}];
        load(fullfile(dirname,'CA1data.mat'),'R','opto');                   % Load rates
        load(fullfile(dirname,'PY_IN.mat'),'PY_IN');                        % Load cell types
        load(fullfile(dirname,'Mcells.mat'),'Mcells');                      % Load mcells
        
        if optnum == 0
            optotrials = (opto(:,2) == optnum);                             % KEEP ALL TRIALS WITHOUT LIGHT 
        else
            optotrials = (opto(:,1) == optnum & opto(:,2) == optoON);       % KEEP ONLY TRIALS OF THAT OPTO WITH LIGHT ON/OFF
        end
        
        R = R(PY_IN == ct, trialbins, optotrials);                          % Keep rates of celltype, across selected bins, and on those trials only                
        
        R = smoothdata(R,2,'gaussian',5);                                   % SMOOTH rates (LESS SMOOTHING THAT OTHER FIGURES!)
        
        % REMOVE TRIALS WHERE IT NEVER SPIKED (TO ELIMINATE DRIFTING UNIT
        % EFFECT) - DOES NOT AFFECT RESULTS
%         for c = 1:size(R,1)
%             nospiketrials = (sum(R(c,:,:),2) == 0);
%             R(c,:,nospiketrials) = nan;
%         end
        
        mr = nanmean(R,3);                                                  % Get mean rates of each cell over all selected trials
        sr = std(R,[],3) / sqrt(size(R,3));                                 % Get SE signal over all trials

        % NORMALIZE BY BASELINE
        mr = mr ./ nanmean(mr(:,baseline),2);                               % Normalize by mean baseline (baseline bins included in trialbins)
        sr = sr ./ nanmean(mr(:,baseline),2);
        mr(isnan(mr) | isinf(mr)) = 0;
        sr(isnan(sr) | isinf(sr)) = 0;
        
        mc = Mcells(ct,:);                                                  % Keep mcells of that celltype
           
        f = nan(size(R,1),2);                                               % Create fields matrix for all cells
        si = nan(size(R,1),1);                                              % and SI matrix
        for tt = 1:2                                                        % For each first odor
            if isempty(mc{tt}), mc{tt} = zeros(0,4);end
            
            f(mc{tt}(:,1),1) = tt;                                          % Tag that odor only in the indexes of mcells
            f(mc{tt}(:,1),2) = mc{tt}(:,2);                                 % and tag the second column with their fields

            si(mc{tt}(:,1)) = mc{tt}(:,4);                                  % Store the SI of Mcells
        end

        if sum(optotrials) == 0
            mr = [];
            sr = [];
            f = [];
            si = [];
        end
        
        MR{ct} = [MR{ct}; mr];                                              % Store mean signal traces
        SR{ct} = [SR{ct}; sr];                                              % Store mean signal traces
        MC{ct} = [MC{ct}; f];                                               % Store fields
        SIn{ct} = [SIn{ct}; si];                                            % Store SI
    end
    
%     zMR{ct} = zscore(MR{ct},[],2);                                           % ZSCORE THE MEAN SIGNALS
    
    % SORT BY (NON ZSCORED) MEAN SIGNALS DURING ODOR
    modor = nanmean(MR{ct}(:,odor1),2);                                     % Mean signal over odor bins
    if isempty(order{ct})
        [~,order{ct}] = sort(modor,'descend');                              % Sort the mean signals over odors (descending order)
    end
    MR{ct} = MR{ct}(order{ct},:);                                           % Reorder accordingly
    SR{ct} = SR{ct}(order{ct},:);                                               
    MC{ct} = MC{ct}(order{ct},:);
    SIn{ct} = SIn{ct}(order{ct});
    % -----------------
    
    zMODOR{ct} = nanmean(MR{ct}(:,odor1),2);                                % Mean zscored signal over the odor
    
    Nr(ct) = size(MR{ct},1);
end

%% COUNT TIME CELLS WITH NEGATIVE ODOR RESPONSES
for ct = 1:2
    z(ct) = 100 * sum(zMODOR{ct} < 1) / length(zMODOR{ct});                 % Compute % of all cells with negative odor response
    
    f = MC{ct}(:,2);                                                        % All fields from all Mcells
    dmc = find(f > 2);                                                      % Find indexes of delay cells
    zz = zMODOR{ct}(dmc);                                                   % Keep their relative odor responses  
    z_de(ct) = 100 * sum(zz < 1) / length(zz);                              % Compute % of time cells with negative odor response
    
    mcells_percentage = sum(~isnan(MC{ct}(:,1)))/size(MC{ct},1)
end
z
z_de

%% PLOT
figure;

for ct = 1:2                                                                % For each cell type
    % PLOT MEAN ZSCORED SIGNALS
    subplot(5,8,(ct-1)*4 + [9:11 17:19 25:27]); hold on;
    maxcol = 2*(ct==1) + 3*(ct==2);
    plot_seq_rates(MR{ct},[],bins,tp,maxcol);
    for tt = 1:2
        k = find(MC{ct}(:,1) == tt);
        plot(MC{ct}(k,2),k,'ok','MarkerFaceColor',C(tt,:),'Markersize',3);
    end
    xlim([0,tp(6)]);
    ylabel(Ctypes{ct});
 
    % PLOT MEAN ZSCORED SIGNALS ACROSS CELLS
    subplot(5,8,(ct-1)*4 + [1:3]); hold on;
    fill_plot(bins,nanmean(MR{ct},1),std(MR{ct},[],1)/sqrt(Nr(ct)),'k');
    axis tight;
    xlim([0,tp(6)]);
    
    % PLOT MEAN ZSCORED SIGNALS ACROSS CELLS WITH POSITIVE/NEGATIVE ODOR RESPONSE
    subplot(5,8,(ct-1)*4 + [1:3]); hold on;
    k = zMODOR{ct}>1;
    fill_plot(bins,nanmean(MR{ct}(k,:),1),std(MR{ct}(k,:),[],1)/sqrt(sum(k)),'b');
    k = zMODOR{ct}<=1;
    fill_plot(bins,nanmean(MR{ct}(k,:),1),std(MR{ct}(k,:),[],1)/sqrt(sum(k)),'r');
    axis tight;
    xlim([0,tp(6)]);

    % PLOT MEAN SIGNALS ACROSS CELLS ZOOMED AROUND ODOR
    subplot(5,8,(ct-1)*4 + 4); hold on;
    k1 = find(bins >= 0.5,1,'first');
    k2 = find(bins <= 2.5,1,'last');
    fill_plot(bins(k1:k2),nanmean(MR{ct}(:,k1:k2),1),std(MR{ct}(:,k1:k2),[],1)/sqrt(Nr(ct)),'k');
    k = zMODOR{ct}>1;
    fill_plot(bins(k1:k2),nanmean(MR{ct}(k,k1:k2),1),std(MR{ct}(k,k1:k2),[],1)/sqrt(sum(k)),'b');
    k = zMODOR{ct}<=1;
    fill_plot(bins(k1:k2),nanmean(MR{ct}(k,k1:k2),1),std(MR{ct}(k,k1:k2),[],1)/sqrt(sum(k)),'r');
    axis tight;
    xlim([bins(k1), bins(k2)]);
    set(gca,'XTick',0.6:0.2:2.4);


    % PLOT MEAN ZSCORED SIGNAL OVER THE ODOR
    subplot(5,8,(ct-1)*4 + [12 20 28]); hold on; view(90, 90);
    plot(1:Nr(ct),zMODOR{ct},'.','Markersize',4);
    line([1 Nr(ct)],[1 1],'LineStyle','--','Color','k');
    axis tight;
    
    if Nr(ct) > 3
        binstep = 3*(ct==1) + 2*(ct==2);
        edges = 1 : binstep : Nr(ct);                                                   % Make bins
        de = edges(2)-edges(1);                                                 % bin size
        [M,S] = bin_x_axis(1:Nr(ct),zMODOR{ct},edges);                              % Bin the MEAN ZSCORED SIGNAL
        fill_plot(edges(1:end-1)+de/2, M, S, 'b');                              % Plot
    end
    
    % PLOT SI
    subplot(5,8,(ct-1)*4 + [33:35]); hold on;
    f = MC{ct}(~isnan(SIn{ct}),2);                                          % Keep all Mcell fields
    si = SIn{ct}(~isnan(SIn{ct}));                                          % and their SI
    plot(f + randn(size(f))*0.06 , si + randn(size(si))*0.01,'.','Markersize',4); % Plot SI as function of field (with jitter)
    
    bon = bins(onbins);                                                     % Keep modulation bins only
    edges = bon(1:20:end);                                                   % Make bins
    de = edges(2)-edges(1);                                                 % bin size
    [M,S] = bin_x_axis(f,si,edges);                                         % Bin the SI distribution
    fill_plot(edges(1:end-1)+de/2, M, S, 'b');                              % Plot
    ylim([0,1.2]); xlim([0.9,7.1]);
    
    % COMPARE SI BETWEEN ODOR ONSET AND REST
    thresh  = 0.2;                                                          % Split at 200ms after odor onset
    if Nr(ct) > 3
        subplot(5,8,(ct-1)*4 + 36); hold on;
        si1 = si(f <= tp(1) + thresh);
        si2 = si(f > tp(1) + thresh);
        p = ranksum(si1,si2,'tail','left');                                  % left-tailed two-sample ttest
        
        dp = table(si1',si2','VariableNames',{'Onset','Post'});
        violinplot(dp);
        hold on;
        plot_significance(p,1,2,si1,si2);
        title(num2str(p));
    end
end