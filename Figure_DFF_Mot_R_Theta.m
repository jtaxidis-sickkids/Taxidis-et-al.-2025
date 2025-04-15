function Figure_DFF_Mot_R_Theta(INclass,cell_flag)

if nargin == 1, cell_flag = 'all'; end

%% LOAD POOLED DATA AND STACK IN SINGLE COLUMN
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'DFF','R','AMP','MC','MOT','bins','time','timepoints');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    DFF{1,5} = [];    DFF{2,3} = [];
    R{1,5} = [];    R{2,3} = [];
    AMP{1,5} = [];    AMP{2,3} = [];
    MC{1,5} = [];   MC{2,3} = [];
    MOT{1,5} = [];   MOT{2,3} = [];
end

% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
% if strcmp(INclass,'SOM')
%     DFF(3,3:7) = DFF(3,4:8);
%     MC(3,3:7) = MC(3,4:8);
% end

DFF = vertcat(DFF{:});
R = vertcat(R{:});
AMP = vertcat(AMP{:});
MC = vertcat(MC{:});
MOT = vertcat(MOT{:});

ls = length(DFF);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
DFFall = {};
Rall = {};
AMPall = {};
MCall = [];
MOTall = {};

for i = 1:ls
    for r = 1:size(DFF{i},1)
        DFFall{count,1} = squeeze(DFF{i}(r,:,:))';                          % Keep only that cell's DFF [trials x time]
        Rall{count,1} = squeeze(R{i}(r,:,:))';                              % Keep only that cell's Rate [trials x bins]
        AMPall{count,1} = squeeze(AMP{i}(r,:,:,2))';                        % Keep only that cell's THETA [trials x time]
        MCall(count,:) = MC{i}(r,:);
        MOTall{count,1} = MOT{i};                                           % Keep only that cell's Locomotion [time x trials]
        count = count + 1;
    end
end
DFF = DFFall;
R = Rall;
AMP = AMPall;
MC = MCall;
MOT = MOTall;
clear DFFall Rall AMPall MCall MOTall

%% CHOOSE CELL GROUP TO USE
Nr = size(MC,1);                            % Number of all cells
if strcmp(cell_flag,'all')
    k = true(Nr,1);
elseif strcmp(cell_flag,'mcells')
    k = ~isnan(MC(:,1));
elseif strcmp(cell_flag,'odorcells')
    k = MC(:,1) > 0;
elseif strcmp(cell_flag,'nonodorcells')
    k = MC(:,1) == 0;
elseif strcmp(cell_flag,'nonmcells')
    k = isnan(MC(:,1));
end

DFF = DFF(k);
R = R(k);
AMP = AMP(k);
MOT = MOT(k);

Nr = sum(k);
Ntr = cellfun(@(x) size(x,1),DFF);

%% SMOOTH RAW DFF
for i = 1:Nr
    DFF{i} = sgolayfilt(DFF{i}, 1, 21, [] , 2);                             % Apply linear filter along time bins
end

%% GET HYPERPOLARIZATIONS + THEIR ONSET + DURATION
minD = cell(Nr,1);
minON = cell(Nr,1);
minDUR = cell(Nr,1);
for c = 1:Nr                                                                % For each selected cell
    [minD{c},~,minON{c},minDUR{c}] = get_DFF_minmax(DFF{c},time,timepoints);% NORMALIZE ODOR DFF TO BASELINE AND COMPUTE MIN/MAX
    minD{c} = abs(minD{c});                                                 % Keep absolute value of hyperpolarization
end

%% SET TIME FOR MOTION/THETA AND BINS FOR RATE COMPUTATION
minT = cell(Nr,1);
minB = cell(Nr,1);
for c = 1:Nr                                                                % For each selected cell
    minT{c} = [minON{c}, minON{c} + minDUR{c}];                             % Compute hyperpolarization onset/offset times (time starts at 0)
    minT{c} = minT{c} + 1;                                                  % Add + 1 to switch to trial time
    
    minB{c} = cell(Ntr(c),1);                                               
    for tr = find(~isnan(minON{c}))'                                        % For all trials with a significant hyperpol.
        minB{c}{tr} = find(bins > minT{c}(tr,2) & bins < 2);                % Find the bins AFTER the hyperpol and before odor offset.
    end
    
    minT{c} = round(minT{c} * 1000);                                        % Turn hyperpol. trial time to indexes on time vector
end

%% COMPUTE HYPERPOL AMPLITUDE AND FREQUENCY OVER ALL TRIALS    
MminD = cellfun(@nanmean, minD);                                            % Get average hyperpolarization (non-significant ones are not accounted for)  
    
FminD = nan(Nr,1);
for c = 1:Nr                                                                % For each cell
    FminD(c) = 100 * sum(~isnan(minD{c})) / length(minD{c});                % Count ratio of significant hyperpolarizations over trials
end

%% COMPUTE LOCOMOTION DURING HYPERPOLARIZATION
hMOT = cell(Nr,1);
for c = 1:Nr                                                                % For each selected cell
    hMOT{c} = nan(Ntr(c),1);                                                % Set vector for all trials
    for tr = find(~isnan(minT{c}(:,1)))'                                    % For all trials with a significant hyperpol.
        hMOT{c}(tr) = mean(MOT{c}(minT{c}(tr,1):minT{c}(tr,2), tr));        % Mean motion during hyperpolarization for each trial
    end
end
mhMOT = cellfun(@nanmean,hMOT);                                             % Average across trials per cell

%% COMPUTE THETA AMPLITUDE DURING HYPERPOLARIZATION
hAMP = cell(Nr,1);
for c = 1:Nr                                                                % For each selected cell
    hAMP{c} = nan(Ntr(c),1);
    for tr = find(~isnan(minT{c}(:,1)))'
        hAMP{c}(tr) = mean(AMP{c}(tr, minT{c}(tr,1):minT{c}(tr,2)));        % Mean theta amplitude during hyperpolarization for each trial
    end
end
mhAMP = cellfun(@nanmean,hAMP);

%% COMPUTE RATES AFTER HYPERPOLARIZATION
hR = cell(Nr,1);
for c = 1:Nr                                                                % For each selected cell
    hR{c} = nan(Ntr(c),1);% Set vector for all trials
    k = cellfun(@isempty, minB{c});                                         % Find trials without significant hyperpol. (OR hyperpol ends at odor offset)
    for tr = find(~k)'                                                      % For all trials with a significant hyperpol.
        hR{c}(tr) = mean(R{c}(tr, minB{c}{tr}));                            % Mean rate AFTER hyperpolarization for each trial
    end
end
mhR = cellfun(@nanmean,hR);


%% PLOT HYPERPOLARIZATION AMPLITUDE VS MEASURES
A = {mhMOT; mhR; mhAMP};
B = {MminD; FminD};
lA = length(A);

figure;
xl = {'Locomotion'; 'Rates'; 'Theta Amp'};
yl = {'Hyperp Ampl'; 'Hyperp Freq'};

for i = 1:lA
    for j = 1
        subplot(lA,2,2*(i-1)+j); hold on;
        scatter(A{i}, B{j},'ok')
        
        % LEAST SQUARES ESTIMATE FOR DISTRIBUTION
        k = ~isnan(A{i});                                                   % REMOVE CELLS WITH NO HYPERPOL.
        coefficients = polyfit(A{i}(k),B{j}(k), 1);                         % Equivalent to: [ones(size(D)), D]\C
        Cfitted = polyval(coefficients, [min(A{i}(k)) max(A{i}(k))]);
        pA = coefTest(fitlm(A{i}(k), B{j}(k),'linear'));
        plot([min(A{i}(k)) max(A{i}(k))] , Cfitted,'-k')
        title(num2str(pA));
        xlabel(xl{i}); ylabel(yl{j});
        axis tight;
        
        mdl = fitlm(A{i},B{j})
    end
end

%% COMPUTE LOCOMOTION, THETA AND RATES DURING HYPERP. TIME WINDOW IN ALL TRIALS
hyperp_time = find(time > 1 & time <= 1.15);                                % Hyperpol time limits for MOTION+DFF+THETA (SAME TIME VECTOR)
hyperp_bins = find(bins > 1.15 & bins <= 2);                                % Post-hyperpol time limits for Rate 

mMOT = cell(Nr,1);
mAMP = cell(Nr,1);
mR = cell(Nr,1);
for c = 1:Nr                                                                % For each selected cell
        mMOT{c} = mean(MOT{c}(hyperp_time, :),1);        % Mean motion during hyperpolarization for each trial
        mAMP{c} = mean(AMP{c}(:, hyperp_time),2);        % Mean theta amplitude during hyperpolarization for each trial
        mR{c}   = mean(R{c}(:, hyperp_bins),2);                            % Mean rate AFTER hyperpolarization for each trial
end
mMOT = cellfun(@nanmean,mMOT);                                             % Average across trials per cell
mAMP = cellfun(@nanmean,mAMP);
mR = cellfun(@nanmean,mR);

%% PLOT HYPERPOLARIZATION FREQUENCY VS MEASURES
A = {mMOT; mR; mAMP};
B = {MminD; FminD};
lA = length(A);

for i = 1:lA
    for j = 2
        subplot(lA,2,2*(i-1)+j); hold on;
        scatter(A{i}, B{j},'ok')
        
        % LEAST SQUARES ESTIMATE FOR DISTRIBUTION
        coefficients = polyfit(A{i},B{j}, 1);                         % Equivalent to: [ones(size(D)), D]\C
        Cfitted = polyval(coefficients, [min(A{i}) max(A{i})]);
        pA = coefTest(fitlm(A{i}, B{j},'linear'));
        plot([min(A{i}) max(A{i})] , Cfitted,'-k')
        title(num2str(pA));
        xlabel(xl{i}); ylabel(yl{j});
        axis tight;
        
        mdl = fitlm(A{i},B{j})
    end
end


%% GET TRIALS WITH HIGHEST AND LOWEST LOCOMOTION DURING HYPEROPOLARIZATION
for c = 1:Nr   
    MOT{c} = smoothdata(MOT{c},1,'movmean',50);
end

allMOT = cell2mat(hMOT);
P = prctile(allMOT,[5 95]);

loMOT = nan(Nr,length(time));
hiMOT = nan(Nr,length(time));
loDFF = nan(Nr,length(time));
hiDFF = nan(Nr,length(time));

for c = 1:Nr   
    k1 = hMOT{c} < P(1);
    k2 = hMOT{c} > P(2);
    
    loMOT(c,:) = nanmean(MOT{c}(:,k1),2)';    
    hiMOT(c,:) = nanmean(MOT{c}(:,k2),2)';
    
    loDFF(c,:) = nanmean(DFF{c}(k1,:),1);
    hiDFF(c,:) = nanmean(DFF{c}(k2,:),1);
end

figure;
for c = 1:Nr 
    subplot(221); hold on
    plot(time,loDFF(c,:),'Color',[0.8 0.8 0.8])
    
    subplot(223); hold on
    plot(time,loMOT(c,:),'Color',[1 0.5 0.5])
    
    subplot(222); hold on
    plot(time,hiDFF(c,:),'Color',[0.8 0.8 0.8])
    
    subplot(224); hold on
    plot(time,hiMOT(c,:),'Color',[1 0.5 0.5])
end

subplot(221);
fill_plot(time,nanmean(loDFF,1),SEM(loDFF),'k');
xlim([0 3]);

subplot(223);
fill_plot(time,nanmean(loMOT,1),SEM(loMOT),'r');
xlim([0 3]);

subplot(222);
fill_plot(time,nanmean(hiDFF,1),SEM(hiDFF),'k');
xlim([0 3]);

subplot(224);
fill_plot(time,nanmean(hiMOT,1),SEM(hiMOT),'r');
xlim([0 3]);
