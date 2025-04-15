function Compare_Opto_Rates(sessions)

%% POOL ODOR FIRING RATES WITH AND WITHOUT OPTO
RON = cell(3,1);
ROFF = cell(3,1);
PYIN = cell(3,1);

for s = 1:length(sessions)
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'R','opto','bins');
    load(fullfile(dirname,'PY_IN.mat'),'PY_IN');
    
    if any(opto(:,1) > 0)
        odorbins = (bins > 1 & bins < 2);
        baseline = bins < 1 ;

        % REMOVE TRIALS WHERE IT NEVER SPIKED (TO ELIMINATE DRIFTING UNIT EFFECT)
        for c = 1:size(R,1)
            nospiketrials = (sum(R(c,:,:),2) == 0);
            R(c,:,nospiketrials) = nan;
        end
        
        Roff = R(:,odorbins,  opto(:,1) == 0);
        RBoff = R(:,baseline, opto(:,1) == 0);
        
        Ron = R(:,odorbins,   opto(:,1) == 4);
        RBon = R(:,baseline,  opto(:,1) == 4);
        
        Roff = nanmean(Roff,[2,3])./nanmean(RBoff,[2,3]);
        Ron = nanmean(Ron,[2,3])./nanmean(RBon,[2,3]);
        
        k = 1*strcmp(sessions{s}(1:2),'CD') + 2*strcmp(sessions{s}(1:2),'ZD');
        
        ROFF{k} = [ROFF{k}; Roff];
        RON{k} = [RON{k}; Ron];
        PYIN{k} = [PYIN{k}; PY_IN];
    end
end

for i =1:2
    k = isinf(RON{i}) | isinf(ROFF{i});
    RON{i}(k) = [];
    ROFF{i}(k) = [];  
    PYIN{i}(k) = [];    
end

RON{3} = cell2mat(RON);
ROFF{3} = cell2mat(ROFF);
PYIN{3} = cell2mat(PYIN);

%% PLOT FOR ALL CELLS OR ONLY PV-CRE OR ONLY SST-CRE
figure;
labels = {'PV-Cre';'SST-Cre';'All'};
for i = 1:3
    subplot(3,2,2*i-1); hold on
    py = PYIN{i} == 1;
    plot([1 2],[ROFF{i}(py)'; RON{i}(py)'],'o-b');
    errorbar([0.9 2.1],[nanmean(ROFF{i}(py)) nanmean(RON{i}(py))],[SEM(ROFF{i}(py)) SEM(RON{i}(py))],'ok-')   % Plot means over all cells
    
    pPY = ranksum(ROFF{i}(py), RON{i}(py)) % ttest(ROFF{i}(py), RON{i}(py))                                            % Compare them
    plot_significance(pPY,0.9,2.1,nanmean(ROFF{i}(py)), nanmean(RON{i}(py)))
    xlim([0.8 2.2]);
    ylabel(labels{i});
    title(num2str(pPY));

    subplot(3,2,2*i); hold on
    in = PYIN{i} == 2;
    plot([1 2],[ROFF{i}(in)'; RON{i}(in)'],'o-r');
    errorbar([0.9 2.1],[nanmean(ROFF{i}(in)) nanmean(RON{i}(in))],[SEM(ROFF{i}(in)) SEM(RON{i}(in))],'ok-')   % Plot means over all cells
    
    pIN = ranksum(ROFF{i}(in), RON{i}(in)) % ttest(ROFF{i}(in), RON{i}(in))                                            % Compare them
    plot_significance(pIN,0.9,2.1,nanmean(ROFF{i}(in)), nanmean(RON{i}(in)))
    xlim([0.8 2.2]);
    title(num2str(pIN));
end

