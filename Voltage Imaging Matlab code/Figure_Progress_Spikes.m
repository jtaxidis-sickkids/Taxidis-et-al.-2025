function Figure_Progress_Spikes(INclass,ton)  

%% LOAD POOLED DATA
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'SP','MC');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    SP{1,5} = [];   SP{2,3} = [];
    MC{1,5} = [];   MC{2,3} = [];
end

% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE 
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
if strcmp(INclass,'SOM')
    SP(3,3:7) = SP(3,4:8);
    MC(3,3:7) = MC(3,4:8);
end

%% LOAD AND ORGANIZE DATA
[Na,maxd] = size(SP);
R = cell(Na,maxd);
MR = cell(Na,maxd);

binlen = 0.005;                                                             % SHORT BIN LENGTH OF 5ms
bins = 0 : binlen : 10;%timepoints(6);
lb = length(bins);

for a = 1:Na
    for d = 1:maxd
        count = 1;
        if ~isempty(SP{a,d})
            for c = 1:size(SP{a,d},1)                                       % For each video
                sp = SP{a,d}{c};                                            % Keep all spikes from all cells
                
                [Nr,Ntr] = size(sp);
                for r = 1:Nr                                                % For each cell in the video
                    R{a,d}{count,1} = nan(Ntr,lb-1);
                    for tr = 1:Ntr                                          % For each trial
                        h = histcounts(sp{r,tr},bins) / binlen;             % COMPUTE FIRING RATE
                        R{a,d}{count,1}(tr,:) = smooth(h,3);                % AND SMOOTH
                    end
                    count = count+1;
                end
            end
        end
    end
end

%% AVERAGE ACROSS TRIALS
for a = 1:Na
    for d = 1:maxd
        if ~isempty(R{a,d})
            MR{a,d} = (cellfun(@(x) mean(x,1), R{a,d}, 'UniformOutput', false)); % Mean rate across all trials for each cell
        end
        MR{a,d} = cell2mat(MR{a,d});                                        % Turn to matrix
    end
end

%% AVERAGE OVER SELECTED BINS
onbins = find(bins >= ton(1) & bins <= ton(2));
    
for a = 1:Na
    for d = 1:maxd
        if ~isempty(R{a,d})
            MR{a,d} = MR{a,d}(:,onbins);                                    % Keep selected bins
            MR{a,d} = mean(MR{a,d},2);                                      % Get mean over bins
        end
    end
end

%% PLOT  ACROSS DAYS AND COMPARE NAIVE-TRAINED
A = cell(1,maxd);
for d = 1:maxd                                                              % For each day
    atemp = MR(:,d);                                                        % Keep all mice
    atemp = cell2mat(atemp);                                            	% Turn to matrix
    A{d} = atemp;                                                           % Store are cell entry
end

k = cellfun(@length,A);                                                     % Count #cells in each day
kmax = max(k);                                                              % Max # cells on a single day
Dd = nan(kmax,maxd);
for d = 1:maxd
    Dd(1:k(d),d) = A{d};
end

[Rs,Ps] = Spearman_over_days(Dd);

figure;
subplot(maxd,4,[2 3 6 7]); hold on;
for d = 1:maxd
    plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , Dd(:,d),'.','Color',[0.8 0.8 0.8]);
end
fill_plot(1:maxd,nanmean(Dd),SEM(Dd),'b')
title([num2str(Rs),' (',num2str(Ps),')']);
ylabel('Mean rate at bins');
xlim([0.5 maxd+0.5]);

subplot(maxd,4, [4 8]);
plot_naive_trained(Dd);
