function Get_Modulation_SpikeSignal(day)

%% GET SESSION FOLDERS
reps = 500;                                                                 % Number of repetition for rates shuffling in each trial

rng(1);                                                                     % For reproducibility
celltypes = {'PY';'IN'};

[~,ses] = fileparts(day);
disp(ses);

%% SPLIT TRIALS BY CRITERION
load(fullfile(day,'Trials_delay_5.mat'),'Trials','time','timepoints','Npy');
load(fullfile(day,'Trials_delay_5_SpikeSignal.mat'),'S');

% Split trials
trials = floor(Trials(:,1)/10);                                         % keep the first odor identity
ontime = (time >= timepoints(1) & time <= timepoints(3));               % Use bins from first odor onset to end of delay
ontime1 = find(ontime,1,'first');                                       % Keep the first bin point to compute modulation over
% --------

lt = sum(ontime);                                                       % Number of bins to be used for modulation calculation

ctypes = 2;                                                             % Number of different cell types (PY and IN)
ttypes = max(trials);                                                   % Number of different trial types
Mcells = cell(ctypes,ttypes);                                           % Allocate memory for storing modulated cells
MR = cell(ctypes,ttypes);                                               % and mean rates of all cells over trial types
SR = cell(ctypes,ttypes);                                               % and rate SD
SI = cell(ctypes,ttypes);                                               % and Selectivity Index

NS = size(S,1);                                                         % Number of ROI (includes PY and IN)
disp(['Total cells = ',num2str(NS)]);

%% KEEP ONLY ONE TYPE OF CELLS
for ct = 1:ctypes                                                       % For each type of cells (PY vs IN)
    Sc = R_ctype(S,ct,Npy);
    Ns = size(Sc,1);
    disp(['Total ',celltypes{ct},' cells = ',num2str(Ns)]);
    
    %% KEEP ONLY TRIALS OF GIVEN TYPE
    for tt = 1:ttypes                                                   % For each type of trials
        ind = (trials == tt);                                           % Find trials of that type
        Stemp = Sc(:,ontime,ind);                                       % KEEP SIGNALS OF THAT CELL TYPE OVER MODULATION BINS ONLY AND ON TRIALS OF THAT TYPE
        mS = mean(Stemp,3);                                             % Get mean signals over all trials along modulation-bins
        [maxS,maxtime] = max(mS,[],2);                                  % Keep maximum mean signal and its bin
        
        %% FIND NUMBER OF ACTIVATED TRIALS OF THAT TYPE
        mS = sum(Stemp,2);                                              % Add signal over all modulation bins for each trial
        mS = squeeze(mS);                                               % Turn to cells x trials
        mS = logical(mS);                                               % Turn to logical for each trial (1 = it spiked at least once, 0 = silent)
        activetrials = sum(mS,2);                                       % Add over trials to get #trials with activation
        
        activetrialsthr = max(sum(ind)*0.1 , 3);                        % Threshold of activated trials = 10% of trials of that type or 3 (whichever is largest)
        
        %% FIND THRESHOLD BY CIRC-SHIFTING EACH TRIAL RATE
        Ntr = sum(ind);                                                 % Number of trials of that type
        maxdist = zeros(Ns,reps);
        for r = 1:reps                                                  % For each repetition
            lags = 2*rand(Ns,Ntr) - 1;                                  % (Ns x Ntr) random numbers in [-1 1] range
            lags = floor(lags*lt/2);                                    % Get random (round) lags up to +- duration/2 of modulation-bins
            
            A = zeros(size(Stemp));
            for c = 1:Ns                                                % For each cell
                for tr = 1:Ntr                                          % For each trial
                    A(c,:,tr) = circshift(Stemp(c,:,tr),lags(c,tr));    % Circularly shift signal (over modulation bins) by lag
                end
            end
            A = mean(A,3);                                              % Compute new mean signal over all trials
            maxdist(:,r) = max(A,[],2);                                 % Get maximum signal over modulation bins for each cell (Ns x 1)
        end
        maxthr = prctile(maxdist,95,2);                                 % 95% percentile of maximum meanr rate for each cell (Ns,1) over all repetitions
        
        %% COMPARE
        mcells = (maxS > maxthr & activetrials >= activetrialsthr)';    % GET CELLS THAT SATISFY BOTH CRITERIA
        mcells = find(mcells);                                          % Find the cell indexes
        maxtime = maxtime(mcells);                                      % Keep max-bins of modulated cells
        
        %% SORT
        [maxtime,order] = sort(maxtime);                                % Sort max-bins to get field times
        mcells = mcells(order)';
        mpeaks = maxS(mcells);                                          % and their max signal
        lmc = length(mcells);
        
        %% COMPUTE SI AND REMOVE NEGATIVE SIs
        indx = ~(trials == tt);                                         % Find trials of the opposite type
        Stemp = Sc(mcells,ontime,indx);                                 % KEEP signal OF THAT CELL TYPE OVER MODULATION BINS ONLY AND ON TRIALS OF THAT TYPE
        mS = mean(Stemp,3);                                             % Get mean rates over modulation-bins
        R2 = zeros(lmc,1);
        for c = 1:lmc                                                   % For each mcell
            R2(c) = mS(c,maxtime(c));                                   % Compute its firing rate at its maxbin in the opposite trials
        end
        si = (mpeaks - R2)./(mpeaks + R2);                              % Compute Selectivity Index
        
        out = si < 0;                                                   % Find cells that prefer the opposite trials
        mcells(out) = [];                                               %
        maxtime(out) = [];                                              % REMOVE THEM
        mpeaks(out) = [];                                               %
        si(out) = [];                                                   %
        
        %% STORE
        mtime  = time(maxtime + (ontime1-1))';                          % Get their actual bin times (counting from first bin in trial)
        
        Mcells{ct,tt} = [mcells, mtime, mpeaks];                        % Sort the modulated cells accordingly and store with fields and max rates
        SI{ct,tt} = si;
        
        MR{ct,tt} = mean(Sc(:,:,ind),3);                                % Store mean Rate of ALL cells of that type over all trials of that type
        SR{ct,tt} = std(Sc(:,:,ind),[],3)/sqrt(Ntr);                    % and standard errors
        
        lmc = length(mcells);                                           % Number of modulated cells
        disp([num2str(lmc),' cells were modulated (',num2str(lmc/Ns*100),'%)']);
    end
end

%% SAVE
modfile = fullfile(day,'Modulation_delay_5_firstodor_ASAP.mat');
save(modfile,'Mcells','MR','SR','SI','trials','ontime','time','timepoints');

%% PLOT SEQUENCES
ctypenames = {'PY';'IN'};

for ct = 1:ctypes
    figure('Name',day);
    for tt = 1:ttypes
        if ~isempty(Mcells{ct,tt})
            mc = Mcells{ct,tt}(:,1);
            title([ctypenames{ct},' firstodor-',num2str(tt)]);
            plot_seq_rates(MR{ct,tt}(mc,:),[],[],[],time,timepoints,2); % Plot mean normalized smooth Rates
            ylabel('Rates')
        end
    end
end
drawnow;

