function Unit_Clustering(sessions)

%% SET UP PARAMETERS
Fs = 30000;                                                                 % Waveforms are given in the original sampling 
time = (-14:25)/Fs * 1000;

Ns = length(sessions);
Nc = zeros(Ns,1);                
Nc_init = zeros(Ns,1);                
OUT = cell(Ns,1);

rng(1);                                                                     % For reproducibility

%% POOL ANIMALS
Peak = [];
Trough = [];
PTdist = [];
HW = [];
mR = [];
mW = [];
BI = [];

for s = 1:Ns
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data_WF_new.mat'),'WF');
    load(fullfile(dirname,'CA1data.mat'),'CA1spikes','R','opto');
    
    peak = WF.Peak;
    trough = WF.Trough;
    ptdist = WF.PTdist;
    hw = WF.HW;
    W = WF.newW;
    
    % KEEP ONLY NO-OPTO TRIALS
    R = R(:,:,~opto(:,2));
    CA1spikes = CA1spikes(:,~opto(:,2));
    
    % REMOVE TRIALS WHERE IT NEVER SPIKED (TO ELIMINATE DRIFTING UNIT EFFECT)
    for c = 1:size(R,1)
        nospiketrials = (sum(R(c,:,:),2) == 0);
        R(c,:,nospiketrials) = nan;
        CA1spikes(c,nospiketrials) = {0};
    end
    
    [Nc_init(s),~,Ntr] = size(R);                                          % number of trials 
    % ---------------------------------
    
    % REMOVE CELLS WITH SPARSE SPIKING
    allspikes = cellfun(@length, CA1spikes);
    Nsp = sum(allspikes,2);
    
    mr = nanmean(R,[2,3]);
    out = isnan(peak) | mr < 0.5; %(Nsp < 5*Ntr); %

    peak(out) = [];
    trough(out) = [];
    ptdist(out) = [];
    hw(out) = [];
    R(out,:,:) = [];
    W(out,:,:) = [];
    CA1spikes(out,:) = [];
     
    OUT{s} = out;
    
    Peak = [Peak; peak];                                                          % Keap all relevant measures
    Trough = [Trough; trough];
    PTdist = [PTdist; ptdist];
    HW = [HW; hw];

    mR = [mR; nanmean(R,[2,3])];
    mW = [mW; nanmean(W,3)];
    
    Nc(s) = size(R,1);
    
    bi = nan(Nc(s),1);
    isi = cellfun(@diff, CA1spikes,'UniformOutput',false);
    for c = 1:Nc(s)
        ISI = cat(1, isi{c,:});
        bursts = [0; ISI < 0.02];                               % Get bursty spikes
        bi(c) = sum(bursts) / Nsp(c) * 100;                 % Burst index = percentage of spikes closer than 20 ms
    end
    BI = [BI; bi];
end
BI = log10(BI);
BI(isinf(BI) | isnan(BI)) = 0;
mR(isnan(mR)) = 0;

PTamp = Peak - Trough;                                         % Peak-Trough amplitude
PTr = (Peak./Trough);

%% SPLIT INTO PY AND IN
Y = [PTdist,PTr,mR];

v1 = 'Peak-Trough distance';
v2 = 'Peak-Trough ratio';
v3 = 'mean Rate';

Yz = zscore(Y,[],1);                                                        % z-score each variable
[~,score,~,~,explained] = pca(Yz);                                          % PCA 
Ypca = score(:,1:3);                                                        % Keep first 3 PCs (ALL OF THEM)
expl = explained(1:3)
[idx,C] = kmeans(Ypca,2,'Replicates',100,'Start','plus');                   % Cluster PCs into 2 cluster with kmeans

[~,smallclust] = min([sum(idx==1),sum(idx==2)]);                            % Find the smallest cluster
bigclust = setdiff(1:2,smallclust);
C1 = find(idx == bigclust);
C2 = find(idx == smallclust);
idx(C1) = 1;                                                                % SET THE LARGER CLUSTER AS CLUSTER 1
idx(C2) = 2;

%% PLOT PCA IN MORE DETAIL
figure;
subplot(2,4,1); hold on;
plot(Y(C1,1),Y(C1,2),'o','Color','k','Markerfacecolor','k');                % Plot the large cluster first
plot(Y(C2,1),Y(C2,2),'o','Color','r','Markerfacecolor','r');
xlabel(v1);
ylabel(v2);
title('PCA & K-means');

subplot(2,4,5); hold on;
fill_plot(time,nanmean(mW(C1,:)),SEM(mW(C1,:)),'k')
plot(time,nanmean(mW(C1,:)),'k','linewidth',2);
fill_plot(time,nanmean(mW(C2,:)),SEM(mW(C2,:)),'r')
plot(time,nanmean(mW(C2,:)),'r','linewidth',2);
xlabel('time (ms)');
ylabel('Average waveforms');
lc1 = length(C1);
lc2 = length(C2);
title(['PYs = ',num2str(lc1),' (',num2str(100*lc1/(lc1+lc2)),'%). INs = ',num2str(lc2),' (',num2str(100*lc2/(lc1+lc2)),'%']);
axis tight

subplot(2,4,2); hold on;
plot3(Y(C1,1),Y(C1,2),Y(C1,3),'ok','markerfacecolor','k');
plot3(Y(C2,1),Y(C2,2),Y(C2,3),'or','markerfacecolor','r'); view(3); axis square; grid on;
title('PCA & K-means'); xlabel(v1);ylabel(v2);zlabel(v3);

subplot(2,4,3);hold on;
plot3(Ypca(C1,1),Ypca(C1,2),Ypca(C1,3),'ok','markerfacecolor','k');
plot3(Ypca(C2,1),Ypca(C2,2),Ypca(C2,3),'or','markerfacecolor','r');  view(3); axis square; grid on;
title('PCA space'); xlabel('PCA1');ylabel('PCA2');  zlabel('PCA3');

subplot(2,4,6); hold on;
compare_hist(PTdist,C1,C2);
xlabel('Peak-Trough distance (ms)'); 

subplot(2,4,7); hold on;
compare_hist(PTr,C1,C2);
xlabel('Peak-Trough Ratio'); axis tight;

subplot(2,4,8); hold on;
compare_hist(HW,C1,C2);
xlabel('Half Width (ms)'); axis tight;

%% SPLIT BACK INTO SESSIONS AND SAVE
idxtemp = idx;                                                              % Reproduce the cluster indexes

for s = 1:Ns
    PY_IN = nan(Nc_init(s),1);
    PY_IN(~OUT{s}) = idxtemp(1:Nc(s));                      % Store the cluster indexes of the first Nu units (equal to the units of that set)
    idxtemp(1:Nc(s)) = [];                                          % Delete the corresponding indexes
    
    dirname = ['../ProcessedData/',sessions{s}];
    save(fullfile(dirname,'PY_IN.mat'),'PY_IN')
end

%% COMPARE FIRING RATES WITHOUT OPTO
RPY = [];
RIN = [];
for s = 1:Ns
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'R','opto');
    load(fullfile(dirname,'PY_IN.mat'),'PY_IN');
 
    R = R(:,:,opto(:,2) == 0);
    
    % REMOVE TRIALS WHERE IT NEVER SPIKED (TO ELIMINATE DRIFTING UNIT EFFECT)
    for c = 1:size(R,1)
        nospiketrials = (sum(R(c,:,:),2) == 0);
        R(c,:,nospiketrials) = nan;
    end
    
    RPY = [RPY; nanmean(R(PY_IN == 1,:,:),[2,3])];
    RIN = [RIN; nanmean(R(PY_IN == 2,:,:),[2,3])];
end

% Plot
subplot(2,4,4);hold on;
mr = table(RPY',RIN','VariableNames',{'PY','IN'});
violinplot(mr);
pmr = ranksum(RPY,RIN);
plot_significance(pmr,1,2,RPY,RIN)
xlim([0.8 2.2]);

return

%% COMPUTE AND PLOT ISI AUTOCORRELOGRAM
edges = -0.2:0.002:0.2;
cols = [0 0.4470 0.7410;
        0.8500 0.3250 0.0980
        0 0 0];
    
for s = 10:11%Ns
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'CA1spikes','opto');
    load(fullfile(dirname,'PY_IN.mat'),'PY_IN');
    PY_IN(isnan(PY_IN)) = 3;
 
    % From sessions with opto keep only initial no-opto trials
    optotrials = opto(:,1) > 0;
    CA1spikes(:,optotrials) = [];
    
    [Nc,Ntr] = size(CA1spikes);
    ISIall = cell(Ntr,1);
    ISIcorr = nan(Nc,length(edges)-1);
    
    for c = 1:Nc
        for tr = 1:Ntr
            A = CA1spikes{c,tr} - CA1spikes{c,tr}';
            A = A(~eye(size(A)));                       % REMOVE DIAGONAL AND MAKE COLUMN
            A(abs(A) > 0.2) = [];
            ISIall{tr} = A;
        end
        isi = cell2mat(ISIall);
        ISIcorr(c,:) = histcounts(isi,edges);
    end
    
    allspikes = cellfun(@length, CA1spikes);
    Nsp = sum(allspikes,2);
  
    figure;
    for c = 1:Nc
        subplot(9,9,c);
        bar(edges(1:end-1), ISIcorr(c,:),'EdgeColor','none','FaceColor',cols(PY_IN(c),:))
        axis tight;
        title(num2str(Nsp(c)));
        drawnow;
    end
end

%%
%[WF.PTdist,PTr,WF.HW] and similarly [mR,PTr,WF.HW] 
%maybe second [WF.PTdist,WF.HW,PTamp]

% are better than 

%[WF.PTdist,PTr,ISI]; 
%[WF.PTdist,ISI,WF.HW]
% [WF.PTdist,WF.HW,mR]
% [WF.PTdist,PTr,mR] (small difference though)










