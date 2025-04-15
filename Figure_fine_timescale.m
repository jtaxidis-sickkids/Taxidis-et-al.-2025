function [MR,bins] = Figure_fine_timescale(INclass,t_firstspike,t_rebound)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
elseif strcmp(INclass,'SOM')
    Days_SOM;
end

ylabels = {'OdorA';'OdorB';'Non-spec';'Non-field'};
tits = {'OdorA','OdorB','All trials'};

%% LOAD POOLED DATA AND STACK IN SINGLE COLUMNS
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'SP','MC','TR','timepoints');

SP = vertcat(SP{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

ls = length(SP); % Number of sessions

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
SPall = {};
MCall = [];
TRall = {};

for i = 1:ls
    for r = 1:size(SP{i},1)
        SPall{count,1} = SP{i}(r,:);
        MCall(count,:) = MC{i}(r,:);
        TRall{count,1} = TR{i};
        count = count + 1;
    end
end
SP = SPall;
MC = MCall;
TR = TRall;
clear SPall MCall TRall

%% SPLIT INTO CELL GROUPS AND TRIAL TYPES
SPpool = cell(4,3);

for tt1 = 1:4                                                               % For each cell group
    k = (MC(:,1) == tt1);
    if tt1 == 3
        k = (MC(:,1) == 0);
    elseif tt1 == 4
        k = isnan(MC(:,1));
    end
    sp = SP(k);                                                             % Keep spikes of cell group
    tr = TR(k);                                                             % And their trials
    
    for i = 1:length(sp)                                                    % For each cell
        for tt2 = 1:2                                                       % for each first odor
            sptt = sp{i}(tr{i} == tt2);                                     % keep corresponding trials
            SPpool{tt1,tt2} = [SPpool{tt1,tt2}; {sptt}];                    % and add them to pooled data
        end
        SPpool{tt1,3} = [SPpool{tt1,3};sp(i)];
    end
end

%% SORT FIELDS
Fp = cell(4,1);
for tt1 = 1:4                                                               % For each cell group
    k = (MC(:,1) == tt1);                                                   % Keep corresponding cells
    if tt1 == 3
        k = (MC(:,1) == 0);
    elseif tt1 == 4
        k = isnan(MC(:,1));
    end
    [Fp{tt1},order] = sort(MC(k,2));                                        % Sort fields
    
    for tt2 = 1:3                                                           % For each trial type
        SPpool{tt1,tt2} = SPpool{tt1,tt2}(order);                           % Sort their spikes
    end
end

%% TURN SPIKETIMES TO FINE-TIMED RATES PER CELL
binlen = 0.005;                                                             % SHORT BIN LENGTH OF 5ms
bins = 0 : binlen : 10;%timepoints(6);
lb = length(bins);

R = cell(4,3);
for tt1 = 1:4                                                               % For each cell group
    Nr = size(SPpool{tt1,1},1);
    for tt2 = 1:3
        R{tt1,tt2} = cell(Nr,1);
        for r = 1:Nr
            Ntr = length(SPpool{tt1,tt2}{r});
            R{tt1,tt2}{r} = nan(Ntr,lb-1);
            for tr = 1:Ntr
                h = histcounts(SPpool{tt1,tt2}{r}{tr},bins) / binlen;       % COMPUTE FIRING RATE
                R{tt1,tt2}{r}(tr,:) = smooth(h,3);                          % AND SMOOTH
            end
        end
    end
end

%% GET MEAN SPIKE HISTOGRAMS
MR = cell(4,3);
for tt1 = 1:4                                                               % For each cell type
    for tt2 = 1:3                                                           % For each trial type
        h = cellfun(@(x) mean(x,1), R{tt1,tt2},'UniformOutput',false);      % GET CELL'S MEAN RATE OVER ALL TRIALS OF THAT TYPE
        MR{tt1,tt2} = cell2mat(h);
    end
end

%% SIGNIFICANCE ACROSS ODOR INTERVALS OVER DIFFERENT ODORS
odbins = find(bins >= timepoints(1) & bins < timepoints(2));
sig = nan(4,length(odbins));
Pod = nan(4,1);

for tt1 = 1:4                                                               % For each cell type
    for b = 1:length(odbins)
        sig(tt1,b) = ranksum(MR{tt1,1}(:,odbins(b)),MR{tt1,2}(:,odbins(b)));      % GET CELL'S MEAN RATE OVER ALL TRIALS OF THAT TYPE
    end
    [~,sig(tt1,:)] = fdr(sig(tt1,:));
    
    % SIGNIFICANCE ACROSS THE ENTIRE ODOR PERIOD
    r1 = MR{tt1,1}(:,odbins);
    r2 = MR{tt1,2}(:,odbins);
    [Pod(tt1),ttype(tt1)] =  significance(r1(:),r2(:),'unequal');
end
[~,Pod] = fdr(Pod)
ttype

%% PLOT ALL MCELLS AND NON-MCELLS OVER ALL TRIALS IN ONE PLOT (FOCUS ON 1ST ODOR)
figure; hold on;
count = 0.;

% POOL MCELLS OVER ALL TRIALS
SPmcells = [SPpool{1,3};SPpool{2,3};SPpool{3,3}];
SPnmcells = SPpool{4,3};

% PLOT MCELLS
Nr = length(SPmcells);
cols = parula(round(Nr*1.2));

for r = 1:length(SPmcells)
    for tr = 1:length(SPmcells{r})
        sp = SPmcells{r}{tr};
        plot(sp,ones(size(sp))*count,'.','Color',cols(r,:),'MarkerFaceColor',cols(r,:),'Markersize',2);
        count = count + 0.2;
    end
    count = count + 0.3;
end

% PLOT NON MCELLS
Nr = length(SPnmcells);
cols = gray(round(Nr*1.2));

for r = 1:Nr
    for tr = 1:length(SPnmcells{r})
        sp = SPnmcells{r}{tr};
        plot(sp,ones(size(sp))*count,'.','Color',cols(r,:),'MarkerFaceColor',cols(r,:),'Markersize',2);
        count = count + 0.2;
    end
    count = count + 0.3;
end

axis tight;
xlim([0.8, 2.2]);
ylabel('MCells & Non-field');
title('All trials');

%% PLOT AVERAGE FIRING RATES OF ODOR-SPECIFIC AND NON-SPECIFIC MCELLS
figure;
cols = cell_colors;
for tt1 = 1:4                                                  % For each trial type
    for tt2 = 1:2                                              % For each trial type (again)
        subplot(5,2,(tt1-1)*2 + 1); hold on;        
        histogram('BinEdges',bins,'BinCounts',mean(MR{tt1,tt2},1),...
            'EdgeColor','none','FaceColor',cols(tt2,:));
    end
    for tt2 = 1:2                                                           % Do a separate loop to plot on top of histograms
        fill_plot(bins(1:end-1),mean(MR{tt1,tt2},1),SEM(MR{tt1,tt2}),cols(tt2,:)); % Plot mean+SEM curves
    end
    plot(bins(odbins(sig(tt1,:) < 0.05)), 30*ones(1,sum(sig(tt1,:)<0.05)),'.k')
    
    axis tight;
    xlim([0.8 2.2]);
    ylabel(ylabels{tt1});
    if tt1 == 1, title(tits{tt2}); end
    
    % PLOT OVER ALL TRIALS
    subplot(5,2,(tt1-1)*2 + 2); hold on;
    histogram('BinEdges',bins,'BinCounts',mean(MR{tt1,3},1),...
        'EdgeColor','none','FaceColor',cols(3,:));
    fill_plot(bins(1:end-1),mean(MR{tt1,3},1),SEM(MR{tt1,3}),cols(3,:));
    axis tight;
    xlim([0.8 2.2]);
    if tt1 == 1, title(tits{3}); end
end

%% PLOT AVERAGE FIRING RATES OF ODOR-SPECIFIC MCELLS FOR PREFERRED VS NON-PREFERRED
MRpref = [MR{1,1}; MR{2,2}];
MRnpref = [MR{1,2}; MR{2,1}];

sig = ones(1,length(odbins));
for b = 1:length(odbins)
    sig(b) = ranksum(MRpref(:,odbins(b)),MRnpref(:,odbins(b)));      % GET CELL'S MEAN RATE OVER ALL TRIALS OF THAT TYPE
end
[~,sig] = fdr(sig);

subplot(5,2,9); hold on;
histogram('BinEdges',bins,'BinCounts',mean(MRpref,1),...
    'EdgeColor','none','FaceColor',[0.8 0.8 0.8]);
histogram('BinEdges',bins,'BinCounts',mean(MRnpref,1),...
    'EdgeColor','none','FaceColor','r');
fill_plot(bins(1:end-1),mean(MRpref,1),SEM(MRpref),[0.8 0.8 0.8]);
fill_plot(bins(1:end-1),mean(MRnpref,1),SEM(MRnpref),'r');
plot(bins(odbins(sig < 0.05)), 30*ones(1,sum(sig <0.05)),'.k')
axis tight;
xlim([0.8 2.2]);
ylabel(ylabels{tt1});
title('Pref vs Non-Pref');

%% SIGNIFICANCE ACROSS THE ENTIRE ODOR PERIOD
r1 = MRpref(:,odbins);
r2 = MRnpref(:,odbins);
[POD,tytpe] =  significance(r1(:),r2(:),'unequal')

%% PLOT AVERAGE FIRING RATES OF POOLED MCELLS
MRmc = cell2mat(MR(1:3,3));                                                 % Pool mean rates of all Mcells over all trials

% PLOT OVER ALL TRIALS
subplot(5,2,10); hold on;
histogram('BinEdges',bins,'BinCounts',mean(MRmc,1),...
    'EdgeColor','none','FaceColor',cols(3,:));
fill_plot(bins(1:end-1),mean(MRmc,1),SEM(MRmc),cols(3,:));
axis tight;
ylabel('All Mcells');
xlim([0.8 2.2]);

%% COMPARE FIRSTSPIKE RATES AND REBOUND RATES OF MCELLS AND NONMCELLS
MRnmc = MR{4,3};                                                            % Keep mean rates of all nonMcells over all trials

figure;
for i = 1:2                                                                 % For either firstspikes or rebounds
    if i == 1
        onbins = find(bins >= t_firstspike(1) & bins <= t_firstspike(2));
    else
        onbins = find(bins >= t_rebound(1) & bins <= t_rebound(2));
    end
        
    MRmc_f = MRmc(:,onbins);                                                % KEep mean rates over selsected bins
    MRnmc_f = MRnmc(:,onbins);
    
    MRmc_f = mean(MRmc_f,2);                                                % Get mean per cell over all the selected bins
    MRnmc_f = mean(MRnmc_f,2);
    
    p = ranksum(MRmc_f,MRnmc_f);
       
    subplot(1,2,i); hold on
    hw = table(MRmc_f',MRnmc_f','VariableNames',{'MC','nonMC'});
    violinplot(hw);
    plot_significance(p,1,2,MRmc_f,MRnmc_f)
end

subplot(121); 
ylabel('Mean rate at first spike');
subplot(122); 
ylabel('Mean rate at rebound');

%% COMPUTE AND PLOT COHERENCE OF EACH CELL ACROSS TRIALS
% CXY = [];
% for tt1 = 3                                                  % For each trial type
%     for tt2 = 3                                              % For each trial type (again)
%         for r = 1:length(R{tt1,tt2})
%             Cxy = 0;
%             Ntr = size(R{tt1,tt2}{r},1);
%             count = 0;
%             for tr1 = 1:Ntr-1
%                 for tr2 = tr1+1:Ntr
%                     hsp = R{tt1,tt2}{r};
%                     [C,~,F] = wcoherence(hsp(tr1,:),hsp(tr2,:),1/binlen);
%                     Cxy = Cxy + C;
%                     count = count + 1;
%                 end
%             end
%             Cxy(isnan(Cxy) | isinf(Cxy)) = 0;
% %             Cxy(Cxy < 0) = 0;
%
%             CXY(:,:,r) = Cxy/count;
%
% %             figure;
% %             imagesc(bins(1:end-1),F,CXY(:,:,r),[0 0.5]);
% %             drawnow
% %             title('Magnitude-Squared Coherence')
% %             xlabel('Frequency (Hz)')
%         end
%     end
% end
%
% figure;
% imagesc(bins(1:end-1),F,mean(CXY,3),[0 0.5]);
% title('Magnitude-Squared Coherence')
% xlabel('Frequency (Hz)')
% disp('done')

%% COMPUTE AND PLOT COHERENCE BETWEEN CELLS
% for tt1 = 3                                                  % For each trial type
%     for tt2 = 3                                              % For each trial type (again)
%         for r1 = 31:35%length(HSP{tt1,tt2})
%             for r2 = 31:35%length(R{tt1,tt2})
%                 Cxy = [];
%                 Ntr1 = size(R{tt1,tt2}{r1},1);
%                 Ntr2 = size(R{tt1,tt2}{r2},1);
%                 for tr1 = 1:Ntr1
%                     for tr2 = 1:Ntr2
%                         [C,~,F] = wcoherence(R{tt1,tt2}{r1}(tr1,:),R{tt1,tt2}{r2}(tr2,:),1/binlen);
%                         Cxy = cat(3,Cxy,C);
%                     end
%                 end
%                 Cxy(isnan(Cxy) | isinf(Cxy)) = 0;
%                 Cxy = mean(Cxy,3);
%                 figure;
%                 imagesc(bins(1:end-1),F,Cxy,[0 0.5]);
%                 drawnow
%                 title('Magnitude-Squared Coherence')
%                 xlabel('Frequency (Hz)')
%             end
%         end
%     end
% end





