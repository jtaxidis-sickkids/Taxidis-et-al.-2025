function [tr,tp] = split_trials(trials,bins,timep,crit)

% db = bins(2)-bins(1);                                                       % Half bin length (50% overlap)
tr = zeros(size(trials,1),1);
combos = trials(:,1);                                                       % The odor combo in each trial
results = trials(:,2);                                                      % The outcome in each trial

switch crit
    case 'all'                                                              % To compute modulation over all trials
        tr = tr + 1;                                                        % All trials equal one
        tp = (bins >= timep(1) & bins <= timep(3)-0.1); % BINS BETWEEN AND INCLUDING THE TWO ODORS                                               % All bins to be used for computing modulation
        
    case 'firstodor'                                                        % To split trials according to first odor
        tr(combos == 11 | combos == 12) = 1;                                % Trials starting with odor A
        tr(combos == 21 | combos == 22) = 2;                                % Trials starting with odor B
        tp = (bins >= timep(1) & bins <= timep(3)-0.1);                   % Use bins from first odor onset to end of delay
                                                                            % (exclude bins starting/ending outside limits)    
    case 'match'
        tr(combos == 11 | combos == 22) = 1;                                % Match trials
        tr(combos == 12 | combos == 21) = 2;                                % Non-match trials 
        tp = (bins >= timep(3) & bins <= timep(5));                      % Use bins from second odor onset to start of response window
                                                                            
    case 'secondodor'                                                       % To split trials according to second odor
        tr(combos == 11 | combos == 21) = 1;                                % Trials with odor A second
        tr(combos == 12 | combos == 22) = 2;                                % Trials with odor B second
        tp = (bins >= timep(3) & bins <= timep(4));                      % Use bins from second odor onset onset to lick onset

    case 'outcomes'                                                         % To split trials according to outcome
        tr = results;                                                       % Trials have same numbering as possible outcomes
        tp = 1:length(bins);                                                % Use all bins
        
    case 'rewards'                                                          % To split trials according to reward
        tr(results == 1) = 1;                                               % Trials with a water reward
        tr(results ~= 1) = 2;                                               % Trials without
        tp = (bins >= timep(5) & bins <= timep(6));                   % Use bins from lick onset to offset
                
    case 'corrects'                                                         % To split trials according to correct/error
        tr(results == 1 | results == 4) = 1;                                % Correct decision trials
        tr(results == 2 | results == 3) = 2;                                % Wrong decision trials
        tp = 1:length(bins);                                  % ???         % Use all bins 

    case 'licks'                                                            % To split trials according to licks
        tr(results == 1 | results == 3) = 1;                                % Trials with licks
        tr(results == 2 | results == 4) = 2;                                % Trials without licks
        tp = (bins >= timep(3) & bins <= timep(6));                   % From second odor offset (lick starts) to offset
    
    case 'combo'
        tr(combos == 11) = 1;                                               % AA trials
        tr(combos == 12) = 2;                                               % AB trials
        tr(combos == 21) = 3;                                               % BA trials
        tr(combos == 22) = 4;                                               % BB trials
        tp = (bins >= timep(3) & bins <= timep(5));   % !!!!                   % Use bins from second odor onset to start of response window
end                                                                         % WHICH EXCLUDES THE EFFECT OF THE ACTUAL REWARD
