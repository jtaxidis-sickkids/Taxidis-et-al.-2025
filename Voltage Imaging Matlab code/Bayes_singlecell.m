function [t_error,tt_accuracy,bins] = Bayes_singlecell(Rm,Rpred,traintrials,predtrials,bins,bout)


%% REMOVE BINS WITH NO SPIKING IN EITHER TRIAL TYPE (UNDECODABLE)
bout = logical(bout(1,:) .* bout(2,:));                                     % bins with no spikes in either trial type
bins(bout) = [];                                                            % Remove them from bins vector ...
Rm(:,bout) = [];                                                            % ...and from the mean firing rates for both odors

bins_vec = [-bins bins];                                                    % Double the bin time (make trial-A negative)
Rm = [Rm(1,:), Rm(2,:)];                                                    % Concatenate mean rates
lb = length(bins);

% SIGMA FOR GAUSSIAN FOR CONTINUITY CONSTRAIN
dt = 0.1;                                                                   % Distance of bin times
lbtot = 60;                                                                 % total number of bins (so that sigma is not affected by removed bins)
sigma = dt*lbtot/2;                                                         % SD for gaussian (= half odor-delay duration)
% -------------------------------

%% COMPUTE PROBABILITY OF TRIAL TYPE
Ptt = zeros(1,2);
for tt = 1:2                                                                % For each trial type
    Ptt(tt) = sum(traintrials == tt)/length(traintrials);                   % Compute ratio of training trials of that type
end
Ptt = [Ptt(1)*ones(1,lb), Ptt(2)*ones(1,lb)];                               % Assign trial type probabilities to the whole bins vector
Ptt = Ptt / sum(Ptt);                                                       % Scale so that it s a probability (total = 1) over all bins

%% DECODING
lpred = length(predtrials);
t_error = nan(lpred,lb);
tt_acc = nan(lpred,lb);

for tr = 1:lpred                                                            % For each trial to decode
    Rpr = Rpred(tr,:);                                                      % keep its firing rates (1 x bins)
    
    tprev = 0;
    
    for b = 1:lb                                                            % for each bin to be decoded
        ni = round(dt*Rpr(b));                                              % Keep its 'spike number' in that bin
        Pred = ((dt*Rm).^ni) / factorial(ni);                               % Compute component of posterior probability        
        Pred = Ptt .* (Pred .* exp(-dt*Rm));                                % Bayesian probability for all bins
        Pred = Pred / sum(Pred);                                            % Normalize to sum to 1 (not necessary)
        
        % APPLY CONTINUITY CONSTRAINT ---------
        if tprev > 0                                                        % If this is not the first bin to be decoded in this trial
            K = exp(-(abs(bins_vec) - abs(bins_vec(tprev))).^2/(2*sigma^2)); % Compute distance gaussian probability
            K = K / sum(K);                                                 % Normalize to sum to 1 (no necessary)
            Pred = Pred .* K;                                               % Mulitply with previous to restrict distance
        end
        % -------------------------------------
        
        Pred(isnan(Pred)) = 0;                                              % Set Nans as zero probability

        [~,t_pred] = max(Pred);                                             % Find the peak location
        tprev = t_pred;                                                     % Store the decoded bin to use in restricting next bin
        
        t_pred = bins_vec(t_pred);                                          % Turn decoded bin to time bin
        
        if t_pred < 0                                                       % If it s negative
            t_pred = -t_pred;                                               % Turn to positive
            tt_pred = 1;                                                    % And set to trial type 1
        else
            tt_pred = 2;                                                    % Else type 2
        end
        
        t_error(tr,b) = bins(b) - t_pred;                                   % Store decoded time error
        tt_acc(tr,b) = (tt_pred == predtrials(tr));                         % And whether trial type was correct
    end
end

%% COMPUTE TRIAL-TYPE ACCURACY PER BIN
tt_accuracy = sum(tt_acc)/lpred*100;
