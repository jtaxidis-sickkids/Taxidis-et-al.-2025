function [Rtrain, Rpred, traintr, predtr] = split_train_pred(R1,R2,trials1,trials2,outcomes1,outcomes2,cutoff1,cutoff2,outcome_flag)

Ntr2 = length(trials2);

%% GET THE TRAINING TRIALS AND FIRING RATES
traintr = trials1(1:cutoff1);                                               % Keep training trials up to cutoff1
traintr = traintr(outcomes1(1:cutoff1) == 1);                               % USE ONLY CORRECT TRIALS

Rtrain = R1(:,:,1:cutoff1);                                                 % Keep rates of training trials
Rtrain = Rtrain(:,:,outcomes1(1:cutoff1) == 1);                             % USE ONLY CORRECT TRIALS

%% GET THE PREDICTION TRIALS AND RATES
if outcome_flag == 1                                                        % If decoding only correct trials
    outc = outcomes2(cutoff2+1 : Ntr2);                                     % Keep outcomes from cutoff2 to end

    predtr = trials2(cutoff2+1 : Ntr2);                                     % Keep prediction trials from cutoff2 to end
    predtr = predtr(outc == outcome_flag);                                  % Keep only correct trials

    Rpred = R2(:,:,cutoff2+1 : Ntr2);                                       % Keep rates of prediction trials from cutoff2 to end
    Rpred = Rpred(:,:,outc == outcome_flag);                                % Keep only correct trials

elseif outcome_flag == -1                                                   % If decoding only error trials
    predtr = trials2;                                                       % Keep all trials
    predtr = predtr(outcomes2 == outcome_flag);                             % and then keep only error ones

    Rpred = R2;                                                             % Keep all trial rates
    Rpred = Rpred(:,:,outcomes2 == outcome_flag);                           % and then keep only error ones

elseif outcome_flag == 0                                                    % If decoding all trials
    predtr = trials2(cutoff2+1 : Ntr2);                                     % Keep prediction trials from cutoff2 to end
    Rpred = R2(:,:,cutoff2+1 : Ntr2);                                       % Keep rates of prediction trials from cutoff2 to end
end

%% IF TRAINING ON HITS AND DECODING CR OR FA TRIALS
if abs(outcome_flag) == 1234
    traintr = trials1(outcomes1 == 1);                                      % USE ONLY CORRECT HITS FOR TRAINING FROM ENTIRE SESSION
    Rtrain = R1(:,:,outcomes1 == 1);                                        % USE ONLY CORRECT HITS FOR TRAINING FROM ENTIRE SESSION
    
    if outcome_flag == 1234                                                 % If decoding only CR trials
        predtr = trials2(outcomes2 == 4);                                   % Keep only correct rejections from entire session
        Rpred = R2(:,:,outcomes2 == 4);                                     % Keep only correct rejections from entire session
        
        lother = sum(outcomes2 == 2);                                       % Count the FAs in the session
        if length(predtr) > lother                                          % If CRs are more than the FAs
            out = randperm(length(predtr),length(predtr)-lother);           % Pick some random CRs
            predtr(out) = [];                                               % And remove them to get equal number with FAs
            Rpred(:,:,out) = [];
        end
        
    elseif outcome_flag == -1234                                            % If decoding only FA trials
        predtr = trials2(outcomes2 == 2);                                   % Keep only false alarms from entire session
        Rpred = R2(:,:,outcomes2 == 2);                                     % Keep only false alarms from entire session
        
        lother = sum(outcomes2 == 4);                                       % Count the CRs in the session
        if length(predtr) > lother                                          % Repeat the random downsampling to get
            out = randperm(length(predtr),length(predtr)-lother);           % equal number of FA and CR in each session
            predtr(out) = []; 
            Rpred(:,:,out) = [];
        end        
    end
end


