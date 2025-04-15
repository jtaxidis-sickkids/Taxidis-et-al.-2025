function Get_Modulation_VIversion(session)

dirname = ['../ProcessedData/',session];
load(fullfile(dirname,'CA1data.mat'),'R','opto','progress','bins');
load(fullfile(dirname,'PY_IN.mat'),'PY_IN');

rng(0); % For shuffling consistency

Mcells = cell(2,1);

%% KEEP ONLY TRIALS WITHOUT OPTO
noopto = ~opto(:,2);

R = R(:,:,noopto);
progress = progress(noopto,:);

%% SMOOTH RATES
% figure;hold on
% plot(bins(1:end-1),R(1,:,10));

R = smoothdata(R,2,'gaussian',12);

% plot(bins(1:end-1),Rs(1,:,10));

%% SPLIT TRIALS BY FIRST ODOR
trials = floor(progress(:,1)/10);

onbins = (bins >= 1 + 0.005 & bins <= 7 - 0.005);     % Keep modulation bins (+/- 1/2 bin length around the edges)
onbins1 = find(onbins,1,'first');                                           % Keep the first bin point to compute modulation over

%% ZSCORE RATES
Ron = R(:,onbins,:);                                                        % KEEP RATES  OVER MODULATION BINS ONLY

zRon = zscore(Ron,[],2);
zRon(isnan(zRon)) = 0;
zRon(isinf(zRon)) = 0;

for ct = 1:2
    zRct = zRon(PY_IN == ct,:,:);
    [Ns,~,~] = size(zRct);
    
    %% SPLIT INTO TWO ODORS AND FIND MAXIMUM-RATES TIME BINS
    R1 = zRct(:,:,trials == 1);                                                 % KEEP RATES OVER TRIALS OF EACH ODOR
    R2 = zRct(:,:,trials == 2);
    
    mR1 = mean(R1,3);                                                           % Get mean rates over all trials along modulation-bins
    mR2 = mean(R2,3);                                                           % Get mean rates over all trials along modulation-bins
    mRall = mean(zRct,3);
    
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
    
    maxbin = zeros(Ns,1);
    maxR = zeros(Ns,1);
    
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
    SI = zeros(Ns,1);
    
    for c = 1:Ns                                                                % For each cell
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
    
    Mcells{ct} = [mcells, maxbin, maxR, SI];                                        % Store the maxBins, maxRates and SIs of ALL CELLS
end
Mcells{1}
Mcells{2}

%% SAVE
save(fullfile(dirname,'Mcells_VI.mat'),'Mcells','trials','onbins');

