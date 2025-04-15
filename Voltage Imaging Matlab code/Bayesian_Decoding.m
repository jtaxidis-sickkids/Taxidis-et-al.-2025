function Bayesian_Decoding(INclass,binflag)

if nargin == 1, binflag = []; end

%% LOAD POOLED DATA
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'R','MC','TR','timepoints','bins','onbins');

% REMOVE MM7-1 5th DAY (MAYBE NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    R{1,5} = [];    %R{2,3} = [];
    MC{1,5} = [];   %MC{2,3} = [];
    TR{1,5} = [];  %TR{2,3} = [];
end

% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
% if strcmp(INclass,'SOM')
%     R(3,3:7) = R(3,4:8);
%     MC(3,3:7) = MC(3,4:8);
%     TR(3,3:7) = TR(3,4:8);
% end

%% STACK IN SINGLE COLUMNS
R = vertcat(R{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

ls = length(R);                                                             % Number of sessions

%% INITIALIZE
rng(1);                                                                     % For reproducibility
reps = 500;

T_error = cell(0,3);
TT_acc = cell(0,3);
Bins = cell(0,3);
T_error_sh = cell(0,3,reps);
TT_acc_sh = cell(0,3,reps);
Bins_sh = cell(0,3,reps);

% USE ONLY ODOR-BINS
if strcmp(binflag,'delaybins')
    onbins = (bins >= timepoints(2) & bins <= timepoints(3));
end

bins = bins(onbins);
lb = length(bins);

count = 1;

%% PERFORM DECODING
for s = 1:ls
    Rs = R{s}(:,onbins,:);                                                  % Keep only selected bins
    trials = TR{s};
    mc = MC{s}(:,1);
    [Nr,~,Ntr] = size(Rs);
    
    % SPLIT TRIALS
    cutoff = floor(2*Ntr/3);                                                % SET THE TRAINING TRIALS TO 2/3 TOTAL
    traintr = trials(1:cutoff);                                             % Keep training trials up to cutoff1
    predtr = trials(cutoff+1 : Ntr);                                        % Keep prediction trials from cutoff2 to end
    
    lpred = length(predtr);
    
    %% BAYESIAN DECODING
    if Ntr >= 10 &&  ...                                                     % At least 10 trials total (so 4 predtrials)
            sum(traintr == 1) >= 2 && sum(traintr == 2) >= 2                % and at least two trials of each odor to train with  
        
        for c = 1:Nr   
            Rsc = squeeze(Rs(c,:,:));                                       % Keep the cell's rates
            Rsc = Rsc';                                                     % [trials x bins]
           
            % FIND TYPE OF CELL
            k = 1*(mc(c) > 0) + 2*(mc(c) == 0) + 3*(isnan(mc(c)));
                        
            %% SMOOTH RATES 
            for tr = 1:Ntr
                Rsc(tr,:) = smooth(Rsc(tr,:),5,'moving');
            end
            
            %% SCALE (CANNOT ZSCORE BECAUSE I NEED POSITIVE RATES TO TURN TO SPIKE NUMBER)
            mm = max(Rsc,[],'all');                                         % Find maximum rate over all trials
            Rsc = Rsc ./ max(Rsc,[],2);                                     % Scale each trial by its max (so maxes at 1)
            Rsc = Rsc*mm;                                                   % Then scale by the total max. So each trial maxes at that number
            Rsc(isnan(Rsc) | isinf(Rsc)) = 0;

            %% SHAPE RATES
            Rtrain = Rsc(1:cutoff,:);                                       % Keep rates of training trials
            Rpred = Rsc(cutoff+1 : Ntr,:);                                  % Keep rates of prediction trials
            
            Rm = zeros(2,lb);
            out = zeros(2,lb);
            for tt = 1:2                                                    % For each trial type
                Rtr = Rtrain(traintr == tt,:);                              % keep only rates at that trial type
                Rm(tt,:) = mean(Rtr,1);                                     % Get their mean
                out(tt,:) = (Rm(tt,:) == 0);                                % Keep bins with no spikes
            end
            
            %% TRAIN SVM ON ODOR CELLS OVER 2/3 TRIALS AND PREDICT REST
            [T_error{count,k},TT_acc{count,k},Bins{count,k}] = Bayes_singlecell(Rm,Rpred,traintr,predtr,bins,out);
            
            %% SHUFFLE DATA
            for r = 1:reps                                                  % For each repetition
                lags = 2*rand(1,lpred) - 1;                                 % (1 x lpred) random numbers in [-1 1] range
                lags = floor(lags*lb/2);                                    % Get random (round) lags up to +- duration/2 of modulation-bins
                
                Rpred_sh = Rpred;
                for tr = 1:lpred                                            % For each prediction trial
                    Rpred_sh(tr,:) = circshift(Rpred_sh(tr,:),lags(tr));    % Random shift prediction trial along time axis
                end
                [T_error_sh{count,k,r},TT_acc_sh{count,k,r},Bins_sh{count,k,r}] = Bayes_singlecell(Rm,Rpred_sh,traintr,predtr,bins,out);
            end
            count = count +1;
        end
    end
end

%% PLOT BAYESIAN PERFORMANCE
figure;
for i = 1:3
    T = T_error(:,i)';                      % Keep only one type of Mcells (in a row)
    TT = TT_acc(:,i)';
    B = Bins(:,i)';
    
    
    Tsh = squeeze(T_error_sh(:,i,:));       % Keep the same type of Mcells
    TTsh = squeeze(TT_acc_sh(:,i,:));
    Bsh = squeeze(Bins_sh(:,i,:));
    
    Tsh = Tsh(:)';                       %  Keep all cells and reps in a row 
    TTsh = TTsh(:)';
    Bsh = Bsh(:)';

    % PLOT TIME ERROR AND TRIAL ACCURACY
    sub1 = subplot(2,3,i);
    sub2 = subplot(2,3,i+3);
    [T,TT,B] = plot_Bayesian(T,TT,B,Tsh,TTsh,Bsh,bins,sub1,sub2);
end

return

%%
% COMPUTE CORRELATIONS
[Rt,Pt] = corr(B',T', 'type','Spearman');                                   % Get the correlation
[Rtt,Ptt] = corr(Btt(~isnan(TT))',TT(~isnan(TT))', 'type','Spearman');      % Get the correlation
subplot(sub1);
title(['Spearman correlation = ',num2str(Rt),' (',num2str(Pt),')']);
subplot(sub2);
title(['Spearman correlation = ',num2str(Rtt),' (',num2str(Ptt),')']);

% Plot mean comparison
subplot(2,3,3);hold on;
data = {T(B <= 2); T(B > 2)};
pval = ones(2,2);
[pval(1,2),testtype_T] = significance(data{1},data{2},'unequal'); testtype_T
plot_mean_SE(data,pval,['b';'k'])

% Plot mean comparison
subplot(2,3,6);hold on;
data = {TT(Btt <= 2); TT(Btt > 2)};
pval = ones(2,2);
[pval(1,2),testtype_TT] = significance(data{1},data{2},'unequal'); testtype_TT
plot_mean_SE(data,pval,['b';'k'])

