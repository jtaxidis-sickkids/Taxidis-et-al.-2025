function SVM_Decoding(INclass,binflag)

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

% USE ONLY NAIVE OR TRAINED (REVISIONS) -------------------------
% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
% if strcmp(INclass,'SOM')
%     R(3,3:7) = R(3,4:8);
%     MC(3,3:7) = MC(3,4:8);
%     TR(3,3:7) = TR(3,4:8);
% end
% 
% R = R(:,1:2);
% MC = MC(:,1:2);
% TR = TR(:,1:2);
% R(:,1:2) = [];
% MC(:,1:2) = [];
% TR(:,1:2) = [];
% --------------------

%% STACK IN SINGLE COLUMNS
R = vertcat(R{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

ls = length(R); % Number of sessions

%% INITIALIZE
rng(1);                     % For reproducibility
reps = 500;

Acc = nan(0,3);
Pacc = nan(0,3);
Acc_sh = nan(0,3,reps);
            
% USE ONLY ODOR-BINS     
if strcmp(binflag,'odorbins')
    onbins = (bins >= timepoints(1) & bins <= timepoints(2));
end

count = 1;

%% PERFORM SVM
for s = 1:ls
    Rs = R{s}(:,onbins,:);                                                  % Keep only selected bins
    trials = TR{s};
    mc = MC{s}(:,1);
    [Nr,~,Ntr] = size(Rs);
        
    % SPLIT TRIALS
    trials(trials == 2) = -1;                                               % Set two classes of trial types (+1,-1)
    cutoff = floor(2*Ntr/3);                                                % SET THE TRAINING TRIALS TO 2/3 TOTAL
    traintr = trials(1:cutoff);                                             % Keep training trials up to cutoff1
    predtr = trials(cutoff+1 : Ntr);                                        % Keep prediction trials from cutoff2 to end
    
    lpred = length(predtr);
    ltrain = length(traintr);
    
    % SVM
    if Ntr >=10 &&  ...                                                     % At least 10 trials total (so 4 predtrials)
            sum(traintr == 1) >= 2 && sum(traintr == -1) >= 2               % and at least two trials of each odor to train with

        for c = 1:Nr
            Rsc = squeeze(Rs(c,:,:));                                       % Keep the cell's rates
            Rsc = Rsc';                                                     % [trials x bins]
           
            % FIND TYPE OF CELL
            k = 1*(mc(c) > 0) + 2*(mc(c) == 0) + 3*(isnan(mc(c)));
                        
            %% SMOOTH RATES 
            for tr = 1:Ntr
                Rsc(tr,:) = smooth(Rsc(tr,:),5,'moving');
            end
            
            %% ZSCORE
            Rsc = zscore(Rsc,[],2);                                         % ZSCORE THE RATE OVER MODULATION BINS ONLY!!!
            Rsc(isnan(Rsc) | isinf(Rsc)) = 0;
            
            %% SPLIT RATES
            Rtrain = Rsc(1:cutoff,:);                                       % Keep rates of training trials
            Rpred = Rsc(cutoff+1 : Ntr,:);                                  % Keep rates of prediction trials from cutoff2 to end
            
            %% TRAIN SVM ON ODOR CELLS OVER 2/3 TRIALS AND PREDICT REST
            [pred,Acc(count,k)] = SVM(Rtrain,Rpred,traintr,predtr);            
            [predtr pred]
                        
            %% SHUFFLE DATA
            for r = 1:reps
                Rt_sh = Rtrain(randperm(ltrain),:);                         % ALSO shuffle training trials
                Rp_sh = Rpred(randperm(lpred),:);                           % Shuffle prediction trials (not enough since is like only shuffling order of predictions (??))
                [~,Acc_sh(count,k,r)] = SVM(Rt_sh,Rp_sh,traintr,predtr);    % THIS CAN BE REPLACED BY JUST COMPUTING THE NEW ERROR
            end
            Pacc(count,k) = significance(squeeze(Acc_sh(count,k,:)),Acc(count,k),'smaller');

            %% SET EMPTY ENTRIES IN NEWLY CREATED ROWS AS NaN (OTHERWISE THEY ARE ZERO AND AFFECT RESULTS)
            Acc(count,setdiff(1:3,k)) = nan;
            Acc_sh(count,setdiff(1:3,k),:) = nan;
            Pacc(count,setdiff(1:3,k)) = nan;
            
            count = count + 1;
        end
    end
end
size(Acc)
%% SIGNIFICANCE AGAINST BASELINE
pb = nan(3,1);
for i = 1:3
    pb(i) = ranksum(Acc(:,i),reshape(Acc_sh(:,i,:),[],1),'alpha',0.05,'tail','right');
end
[~,pb] = fdr(pb)


%% SIGNIFICANCE IN SVM DIFFERENCES
p = nan(3,1);
p(1) = ranksum(Acc(:,1),Acc(:,2));
p(2) = ranksum(Acc(:,1),Acc(:,3));
p(3) = ranksum(Acc(:,2),Acc(:,3));
[~,p] = fdr(p)

%% PLOT
figure;
subplot(211)
A = table(Acc(:,1),Acc(:,2), Acc(:,3),'VariableNames',{'Mcells','Non-spec','Non-Mcells'});
violinplot(A);
hold on;
for i = 1:3
    line(i+[-0.2 0.2], nanmean(Acc_sh(:,i,:),'all')*[1 1],'Color','r');
end
plot_significance(p(1),1,2,Acc(:,1),Acc(:,2))
plot_significance(p(2),1,3,Acc(:,1),Acc(:,3))
plot_significance(p(3),2,3,Acc(:,2),Acc(:,3))

subplot(212)
P = table(Pacc(:,1),Pacc(:,2), Pacc(:,3),'VariableNames',{'Mcells','Non-spec','Non-Mcells'});
violinplot(P);

%% REPEAT AFTER POOLING ALL CELLS (REVISIONS)
Acc_all = Acc(:);
Acc_sh_all = Acc_sh(:);

% SIGNIFICANCE AGAINST BASELINE
pb_all = ranksum(Acc_all,Acc_sh_all,'alpha',0.05,'tail','right')

% SIGNIFICANCE IN SVM DIFFERENCES
p = nan(3,1);
p(1) = ranksum(Acc_all,Acc(:,1));
p(2) = ranksum(Acc_all,Acc(:,2));
p(3) = ranksum(Acc_all,Acc(:,3));
[~,p] = fdr(p)

% PLOT
figure;
subplot(211)
A = table(Acc_all',Acc(:,1)',Acc(:,2)', Acc(:,3)','VariableNames',{'All cells','Mcells','Non-spec','Non-Mcells'});
violinplot(A);
hold on;
line(1+[-0.2 0.2], nanmean(Acc_sh_all,'all')*[1 1],'Color','r');

plot_significance(p(1),1,2,Acc_all,Acc(:,1))
plot_significance(p(2),1,3,Acc_all,Acc(:,2))
plot_significance(p(3),2,3,Acc_all,Acc(:,3))

