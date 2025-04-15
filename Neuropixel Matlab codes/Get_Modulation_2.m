function Get_Modulation_2(session)

dirname = ['../ProcessedData/',session];
load(fullfile(dirname,'CA1data.mat'),'CA1spikes','opto','progress');
load(fullfile(dirname,'PY_IN.mat'),'PY_IN');

rng(0); % For shuffling consistency

Mcells = cell(2,2);

%% NEW BINNED RATES
binlen = 0.1;                                                               % 100msec rate bins
bins = 0:binlen:11;
[Nr,Ntr] = size(CA1spikes);
R = nan(Nr,length(bins)-1,Ntr);

for tr = 1:Ntr
    for c = 1:Nr
        R(c,:,tr) = histcounts(CA1spikes{c,tr},bins);
    end
end
R = R/binlen;

%% KEEP ONLY TRIALS OF NO OPTO
noopto = ~opto(:,2);
R = R(:,:,noopto);
progress = progress(noopto,:);

%% SMOOTH RATES
% figure;hold on
% plot(bins(1:end-1),R(1,:,10));

% R = smoothdata(R,2,'gaussian',30);
% R = smoothdata(R,2,'movmean',5);
% plot(bins(1:end-1),R(1,:,10));
% return

%% SPLIT TRIALS BY FIRST ODOR
trials = floor(progress(:,1)/10);

%% KEEP MODULATION BINS ONLY
onbins = (bins >= 1 + 0.05 & bins <= 7 - 0.05);                           % Keep modulation bins (+/- 1/2 bin length around the edges)
onbins1 = find(onbins,1,'first');                                           % Keep the first bin point to compute modulation over
Ron = R(:,onbins,:);                                                        % KEEP RATES  OVER MODULATION BINS ONLY

%%
for ct = 1:2
    Rct = Ron(PY_IN == ct,:,:);
    
    %% ZSCORE RATES
    zRct = zscore(Rct,[],2);
    zRct(isnan(zRct) | isinf(zRct)) = 0;
    
    %% SPLIT INTO TWO ODORS AND FIND MAXIMUM-RATES TIME BINS
    for tt = 1:2
        ind = (trials == tt);
        zR1 = zRct(:,:,ind);                                                % KEEP RATES OVER TRIALS OF EACH ODOR     
        mR = nanmean(zR1,3);                                                % Get mean rates over all trials along modulation-bins

        [maxR,maxbin] = max(mR,[],2);                                       % Keep maximum mean rate and its bin
        
        %% FIND THRESHOLD BY CIRC-SHIFTING EACH TRIAL RATE
        maxthr1 = random_circshift(zR1);
        
        %% FIND NUMBER OF ACTIVATED TRIALS OF THAT TYPE
        mR = sum(zR1,2);                                                    % Add firing rate over all modulation bins for each trial
        mR = squeeze(mR);                                                   % Turn to cells x trials
        mR = logical(mR);                                                   % Turn to logical for each trial (1 = it spiked at least once, 0 = silent)
        activetrials = sum(mR,2);                                           % Add over trials to get #trials with activation
        
        activetrialsthr = max(sum(ind)*0.1 , 3);                            % Threshold of activated trials = 10% of trials of that type or 3 (whichever is largest)
        
        %% COMPARE AND SET MCELLS
        mcells = (maxR > maxthr1 & activetrials >= activetrialsthr);        % GET CELLS THAT SATISFY CRITERIA
        mcells = find(mcells);                                              % Find the cell indexes
        maxbin = maxbin(mcells);                                            % and their fields
        maxR = maxR(mcells);                                                % (used just to reduce the vector length to correct size)      
        
        %% COMPUTE NON-ZSCORED FIELD RATE AND SI
        Ns = length(mcells);
        SI = zeros(Ns,1);
        for c = 1:Ns                                                        % For each cell
            maxR(c) = nanmean(Rct(mcells(c),maxbin(c),ind));                % Replace with the non-zscored mean rate at field
        
            mR1 = maxR(c);
            mR2 = nanmean(Rct(mcells(c),maxbin(c),~ind));                   % Same for opposite odor trials
            
            SI(c) = (mR1 - mR2)./(mR1 + mR2);                               %  Selectivity Index
        end
        
        %% REMOVE NEGATIVE SI (ONLY FOR PYs)
%         if ct == 1
%             out = SI < 0;                                                       % Find cells that prefer the opposite trials
%             mcells(out) = [];                                                   %
%             maxbin(out) = [];                                                   % REMOVE THEM
%             maxR(out) = [];                                                     %
%             SI(out) = [];                                                       %
%         end
        
        %% ORGANIZE OUTPUT
        maxbin = bins(maxbin + (onbins1-1))';                               % Get their actual bin times (counting from first bin in trial)

        Mcells{ct,tt} = [mcells, maxbin, maxR, SI];                         % Store the maxBins, maxRates and SIs of ALL CELLS
    end
end

celldisp(Mcells')

%% SAVE
save(fullfile(dirname,'Mcells_2.mat'),'Mcells','trials','onbins');

