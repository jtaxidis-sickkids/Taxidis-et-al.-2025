function Compare_Opto(sessions)

%% COMPARE ODOR FIRING RATES WITH AND WITHOUT OPTO
RON = [];
ROFF = [];
PYIN = [];

for s = 1:length(sessions)
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'R','opto','bins');
    load(fullfile(dirname,'PY_IN.mat'),'PY_IN');
    
    if any(opto(:,1) > 0)
        odorbins = (bins >= 1 & bins <= 2);
        baseline = bins < 1 ;

        % REMOVE TRIALS WHERE IT NEVER SPIKED (TO ELIMINATE DRIFTING UNIT EFFECT)
        for c = 1:size(R,1)
            nospiketrials = (sum(R(c,:,:),2) == 0);
            R(c,:,nospiketrials) = nan;
        end
        
        Roff = R(:,odorbins,  opto(:,1) == 0);
        Ron = R(:,odorbins,   opto(:,1) == 4);
        RBoff = R(:,baseline, opto(:,1) == 0);
        RBon = R(:,baseline,  opto(:,1) == 4);
        
        ROFF = [ROFF; nanmean(Roff,[2,3]) ./ nanmean(RBoff,[2,3])];
        RON = [RON; nanmean(Ron,[2,3]) ./ nanmean(RBon,[2,3])];
        PYIN = [PYIN; PY_IN];
    end
end

RON(isinf(RON)) = nan;
ROFF(isinf(ROFF)) = nan;

% Plot separately for PY and IN
figure; 
subplot(121); hold on
py = PYIN == 1;
plot([1 2],[ROFF(py)'; RON(py)'],'o-b');

pPY = ranksum(ROFF(py), RON(py)) % ttest(ROFF(py), RON(py))                                            % Compare them
errorbar([0.9 2.1],[nanmean(ROFF(py)) nanmean(RON(py))],[SEM(ROFF(py)) SEM(RON(py))],'ok-')   % Plot means over all cells
plot_significance(pPY,0.9,2.1,nanmean(ROFF(py)), nanmean(RON(py)))
xlim([0.8 2.2]);

subplot(122); hold on
in = PYIN == 2;
plot([1 2],[ROFF(in)'; RON(in)'],'o-r');

pIN = ranksum(ROFF(in), RON(in)) % ttest(ROFF(in), RON(in))                                            % Compare them
errorbar([0.9 2.1],[nanmean(ROFF(in)) nanmean(RON(in))],[SEM(ROFF(in)) SEM(RON(in))],'ok-')   % Plot means over all cells
plot_significance(pIN,0.9,2.1,nanmean(ROFF(in)), nanmean(RON(in)))
xlim([0.8 2.2]);
