function [MinD,MaxD,MinOn,MinDur,MinT,MinAll] = get_DFF_minmax(dff,time,timepoints)

%% DETECTION THRESHOLDS
crit1 = -1;
crit2 = -3;
dur1 = 50;
dur2 = 10;
distcrit = 20;
durhyp = 200;

%% SHAPE DFF
basedff = dff(:,time >= 0.5 & time <= timepoints(1));                       % BASELINE
mu = mean(basedff,2);                                                       % MEAN ACROSS TIME FOR EACH TRIAL
sigma = std(basedff,[],2);                                                  % STD ACROSS TIME FOR EACH TRIAL
dff = (dff - mu) ./ sigma;                                                  % Z-SCORE

oddff = dff(:,time >= timepoints(1) & time <= timepoints(2));               % KEEP ONLY ODOR BINS

%% INITIALIZE
Ntr = size(oddff,1);
MinD = nan(Ntr,1);
MaxD = nan(Ntr,1);
MinOn = nan(Ntr,1);
MinDur = nan(Ntr,1);
MinT = nan(Ntr,1);

%% COMPUTE FEATURES
for tr = 1:Ntr                                                              % For each trial
    tdip_low = find(oddff(tr,:) < crit1);                                   % find any zscored odor DFFs less than -2   
    % figure;plot(DFF{i}(tr,:));hold on;plot(dff(tr,:));
    % plot(tdip_low+1000,0.5*ones(size(tdip_low)),'.-');
    
    k = find(diff(tdip_low) > distcrit,1,'first');                          % find points corresponding to continuous segments
    if isempty(k), k = length(tdip_low);end
    tdip_low = tdip_low(1:k);                                               % keep only those   
    % plot(tdip_low+1000,zeros(size(tdip_low)),'.-');
    
    if ~isempty(tdip_low) & tdip_low(end)-tdip_low(1) > dur1                % If that segment is at least 50msec long
        tdip_high = find(oddff(tr,tdip_low) < crit2);                       % Find if any dff sub-segments lower than 4
        tdip_high = tdip_high + tdip_low(1) - 1;                            % keep the actual indexes of those segments     
        % plot(tdip_high+1000,-0.5*ones(size(tdip_high)),'.-');
        
        if ~isempty(tdip_high) & tdip_high(end)-tdip_high(1) > dur2           % If that segment is at least 10msec long
            % plot(tdip_high+1000,-ones(size(tdip_high)),'.-');
            
            [MinD(tr),tmin] = min(oddff(tr,tdip_high));                     % Keep the value and timing of the minimum
            MinT(tr) = tmin/1000 + tdip_high(1)/1000;                       % turn it to actual timepoint
            
            MinDur(tr) = (tdip_low(end) - tdip_low(1)) / 1000;              % Keep duration of lower-thresh segment (in sec)
            MinOn(tr) = tdip_low(1)/1000;                                   % Keep initial point of lower-thresh segment
            
            [MaxD(tr)] = max(oddff(tr,tdip_low(end):end));                  % Keep the value of the maximum
        end
    end
end

%% COMPUTE ALL MINIMUM VALUES IRRESPECTIVE OF SIGNIFICANT HYPERPOLARIZATION
MinAll = min(oddff(:,1:durhyp),[],2);                                       % Get minimum value over 1st 200msec in all trials


