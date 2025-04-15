function Get_Modulation(ID,day,session)

rng(0); % For shuffling consistency

%% LOAD DATA
asapfile = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
load(Datafile,'V','B');
disp(asapfile);

pars = B.pars;
trials = B.trials;
R = V.R;
bins = V.bins;

clear V B

%% SPLIT TRIALS BY FIRST ODOR
timepoints = [pars.stim1_on, pars.delay_on, pars.stim2_on, pars.stim2_off, ...
              pars.lick_on, pars.vacuum_on,  pars.trial_dur];

onbins = (bins >= timepoints(1) + 0.05 & bins <= timepoints(3) - 0.05);     % Keep modulation bins (+/- 1/2 bin length around the edges)
onbins1 = find(onbins,1,'first');                                           % Keep the first bin point to compute modulation over

[Nr,~,Ntr] = size(R);

%% SMOOTH RATES
for c = 1:Nr
    for tr = 1:Ntr
        R(c,:,tr) = smooth(R(c,:,tr),5,'moving');
    end
end

%% ZSCORE RATES
Ron = R(:,onbins,:);                                                        % KEEP RATES  OVER MODULATION BINS ONLY 

zRon = zscore(Ron,[],2);
zRon(isnan(zRon)) = 0;
zRon(isinf(zRon)) = 0;

%% SPLIT INTO TWO ODORS AND FIND MAXIMUM-RATES TIME BINS
R1 = zRon(:,:,trials == 1);                                                 % KEEP RATES OVER TRIALS OF EACH ODOR
R2 = zRon(:,:,trials == 2);

mR1 = mean(R1,3);                                                           % Get mean rates over all trials along modulation-bins
mR2 = mean(R2,3);                                                           % Get mean rates over all trials along modulation-bins
mRall = mean(zRon,3);

[maxR1,maxbin1] = max(mR1,[],2);                                            % Keep maximum mean rate and its bin
[maxR2,maxbin2] = max(mR2,[],2);                                            % Keep maximum mean rate and its bin

%% FIND THRESHOLD BY CIRC-SHIFTING EACH TRIAL RATE
maxthr1 = random_circshift(R1);
maxthr2 = random_circshift(R2);

%% COMPARE AND SET MCELLS
mcells1 = (maxR1 > maxthr1);                                                % GET CELLS THAT SATISFY SHUFFLE PEAK CRITERION
mcells2 = (maxR2 > maxthr2);                                                % GET CELLS THAT SATISFY SHUFFLE PEAK CRITERION

mcells = 1*mcells1 + 2*mcells2;                                             % Find the preferred odor of each cell
mcells(mcells == 0) = nan;                                                  % Set non-Mcells (= 0) to Nan
mcells(mcells == 3) = 0;                                                    % Set non-odor-specific Mcells (= 3) to 0
                                                                            % (THIS ASSUMES SAME FIELD TIMEBIN FOR EACH ODOR!!!)
                                                                            
%% SET FIELD AND MAX RATE
mb = [maxbin1 maxbin2];
mr = [maxR1 maxR2];

maxbin = zeros(Nr,1);
maxR = zeros(Nr,1);

for c = find(mcells > 0)'                                                   % For all the odor-specific Mcells
    maxR(c) = mr(mcells(c));                                                % Keep Max Rate at preferred odor
    maxbin(c) = mb(mcells(c));                                              % and field time bin
end

for c = find(isnan(mcells))'                                                % For the non-Mcells
    [maxR(c),maxbin(c)] = max(mRall(c,:));                                  % Keep their Max Rate and bin over all trials
end

for c = find(mcells == 0)'                                                  % For non-odor-specific Mcells
    if abs(diff(mb(c,:))) <= 10                                             % If the peaks are closer than 1 sec
        maxbin(c) = round(mean(mb(c,:)));                                   % Keep the mean time bin between the two fields
        maxR(c) = mRall(c,maxbin(c));                                       % And the Max Rate at that bin over all trials
    else                                                                    % If the peaks are further than 1 sec (2 different fields)
        disp('TWO FIELDS!!');
        [maxR(c),k] = max(mr(c,:));                                         % Keep the highest Max Rat eas the field
        maxbin(c) = mb(c,k);                                                % and the corresponding time bin
    end
end

%% COMPUTE SI
SI = zeros(Nr,1);

for c = 1:Nr                                                                % For each cell
    mR1 = squeeze(R1(c,maxbin(c),:));  
    mR2 = squeeze(R2(c,maxbin(c),:));                            
    mR1 = mean(mR1);
    mR2 = mean(mR2);
    
    SI(c) = (mR1 - mR2)./(mR1 + mR2);                                       %  Selectivity Index 
end

%% SET LOW SI AS NON-ODOR-SPECIFIC
mcells(abs(SI) <= 0.42 & mcells > 0) = 0;

%% ORGANIZE OUTPUT 
maxbin  = bins(maxbin + (onbins1-1))';                                      % Get their actual bin times (counting from first bin in trial)

Mcells = [mcells, maxbin, maxR, SI];                                        % Store the maxBins, maxRates and SIs of ALL CELLS

%% SAVE
modfile = fullfile(path,[videoname,'_Mcells.mat']);
save(modfile,'Mcells','trials','onbins','bins','timepoints');

function maxthr = random_circshift(R)

reps = 1000;
[Nr,lb,Ntr] = size(R);

maxdist = zeros(Nr,reps);

for r = 1:reps                                                  % For each repetition
    lags = 2*rand(Nr,Ntr) - 1;                                  % (Ns x Ntr) random numbers in [-1 1] range
    lags = floor(lags*lb/2);                                    % Get random (round) lags up to +- duration/2 of modulation-bins
    A = zeros(size(R));
    for c = 1:Nr                                                % For each cell
        for tr = 1:Ntr                                          % For each trial
            A(c,:,tr) = circshift(R(c,:,tr),lags(c,tr));        % Circularly shift Rate (over modulation bins) by lag
        end
    end 
    A = mean(A,3);                                              % Compute new mean Rate over all trials
    maxdist(:,r) = max(A,[],2);                                 % Get maximum rate over modulation bins for each cell (Ns x 1)
end
maxthr = prctile(maxdist,95,2);                                 % 95% percentile of maximum meanr rate for each cell (Ns,1) over all repetitions
