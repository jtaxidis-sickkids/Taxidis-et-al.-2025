function Get_Modulation_2Pversion(session)

dirname = ['../ProcessedData/',session];
load(fullfile(dirname,'CA1data.mat'),'R','opto','progress','bins');
load(fullfile(dirname,'PY_IN.mat'),'PY_IN');

rng(0); % For shuffling consistency

Mcells = cell(2,2);

%% KEEP ONLY TRIALS WITHOUT OPTO
noopto = ~opto(:,2);

R = R(:,:,noopto);
progress = progress(noopto,:);

%% SMOOTH RATES
R = smoothdata(R,2,'gaussian',12);

%% SPLIT TRIALS BY FIRST ODOR
trials = floor(progress(:,1)/10);

onbins = (bins >= 1 + 0.005 & bins <= 7 - 0.005);     % Keep modulation bins (+/- 1/2 bin length around the edges)
onbins1 = find(onbins,1,'first');                                           % Keep the first bin point to compute modulation over

%% KEEP ONLY TRIALS OF GIVEN TYPE
for ct = 1:2                                                                % For each cell type
    Rct = R(PY_IN == ct,:,:);                                               % Keep only cells of that type
    [Ns,~,~] = size(Rct);
    for tt = 1:2                                                            % For each first odor
        ind = (trials == tt);                                               % Find trials of that type
        Rtemp = Rct(:,onbins,ind);                                          % KEEP RATES OF THAT CELL TYPE OVER MODULATION BINS ONLY AND ON TRIALS OF THAT TYPE
       
        mR = mean(Rtemp,3);                                                 % Get mean rates over all trials along modulation-bins
        [maxR,maxbin] = max(mR,[],2);                                       % Keep maximum mean rate and its bin
        
        %% FIND NUMBER OF ACTIVATED TRIALS OF THAT TYPE
        mR = sum(Rtemp,2);                                                  % Add firing rate over all modulation bins for each trial
        mR = squeeze(mR);                                                   % Turn to cells x trials
        mR = logical(mR);                                                   % Turn to logical for each trial (1 = it spiked at least once, 0 = silent)
        activetrials = sum(mR,2);                                           % Add over trials to get #trials with activation
        
        activetrialsthr = max(sum(ind)*0.1 , 3);                            % Threshold of activated trials = 10% of trials of that type or 3 (whichever is largest)
        
        %% FIND THRESHOLD BY CIRC-SHIFTING EACH TRIAL RATE
        maxthr = random_circshift(Rtemp);
        
        %% COMPARE
        mcells = (maxR > maxthr & activetrials >= activetrialsthr)';        % GET CELLS THAT SATISFY ALL CRITERIA
        mcells = find(mcells);                                              % Find the cell indexes
        maxbin = maxbin(mcells);                                            % Keep max-bins of modulated cells
        
        %% SORT
        [maxbin,order] = sort(maxbin);                                      % Sort max-bins to get field times
        mcells = mcells(order)';
        mpeaks = maxR(mcells);                                              % and their max-Rate
        lmc = length(mcells);
        
        %% COMPUTE SI AND REMOVE NEGATIVE SIs
        indx = ~(trials == tt);                                             % Find trials of the opposite type
        Rtemp = Rct(mcells,onbins,indx);                                    % KEEP RATES OF THAT CELL TYPE OVER MODULATION BINS ONLY AND ON TRIALS OF THAT TYPE
        mR = mean(Rtemp,3);                                                 % Get mean rates over modulation-bins
        R2 = zeros(lmc,1);
        for c = 1:lmc                                                       % For each mcell
            R2(c) = mR(c,maxbin(c));                                        % Compute its firing rate at its maxbin in the opposite trials
        end
        si = (mpeaks - R2)./(mpeaks + R2);                                  % Compute Selectivity Index
        
%         out = si < 0;                                                       % Find cells that prefer the opposite trials
%         mcells(out) = [];                                                   %
%         maxbin(out) = [];                                                   % REMOVE THEM
%         mpeaks(out) = [];                                                   %
%         si(out) = [];                                                       %
        
        %% STORE
        mbins  = bins(maxbin + (onbins1-1))';                          % Get their actual bin times (counting from first bin in trial)
        
%         k = find(PY_IN ==ct);                                   % Keep indexes of all cells of selected celltype
%         mcells_fullindex = k(mcells);                           % Turn mcells in overall index
        Mcells{ct,tt} = [mcells, mbins, mpeaks, si];                        % Sort the modulated cells accordingly and store with fields and max rates
%         SI{ct,tt} = si;
        
        lmc = length(mcells);                                           % Number of modulated cells
        disp([num2str(lmc),' cells were modulated (',num2str(lmc/Ns*100),'%)']);
    end
end

celldisp(Mcells')

%% SAVE
save(fullfile(dirname,'Mcells.mat'),'Mcells','trials','onbins');

return
%% PLOT SEQUENCES
for ct = 1:2
figure;%('Name',day);
    for tt = 1:2
        subplot(1,2,tt);
        if ~isempty(Mcells{ct,tt})
            mc = Mcells{ct,tt}(:,1);
            %         title([ctypenames{ct},' ',crit,'-',num2str(tt)]);
            MR = mean(R(mc,:,trials == tt),3);
            MR = MR ./ max(MR(:,onbins),[],2);
            plot_seq_rates(MR,[],bins,[1 2 7 8 9 11 11]); % Plot mean normalized smooth Rates
            ylabel('Rates')
            
        end
    end
    drawnow;
end


