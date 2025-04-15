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
    A = nanmean(A,3);                                              % Compute new mean Rate over all trials
%             A = smoothdata(A,2,'movmean',21);

    maxdist(:,r) = max(A,[],2);                                 % Get maximum rate over modulation bins for each cell (Ns x 1)
end
maxthr = prctile(maxdist,95,2);                                 % 95% percentile of maximum meanr rate for each cell (Ns,1) over all repetitions
