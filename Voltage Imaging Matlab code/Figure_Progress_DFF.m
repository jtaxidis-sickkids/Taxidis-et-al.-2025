function Figure_Progress_DFF(INclass,cell_flag)  

if nargin == 1, cell_flag = 'all'; end  

%% LOAD POOLED DATA
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'DFF','MC','time','timepoints');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    DFF{1,5} = [];    DFF{2,3} = [];
    MC{1,5} = [];   MC{2,3} = [];
end

% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE 
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
if strcmp(INclass,'SOM')
    DFF(3,3:7) = DFF(3,4:8);
    MC(3,3:7) = MC(3,4:8);
end

[Na,maxd] = size(DFF);
MD = cell(Na,maxd);
MinD = cell(Na,maxd);
MaxD = cell(Na,maxd);
FminD = cell(Na,maxd);

%% LOAD AND ORGANIZE DATA
for a = 1:Na
    for d = 1:maxd
        if ~isempty(DFF{a,d})
            mcells = cell2mat(MC{a,d});
            if isempty(mcells), mcells = zeros(0,4); end
            Nr = size(mcells,1);                                       % Number of all cells
            
            % CHOOSE CELL GROUP TO USE
            if strcmp(cell_flag,'all')
                k = true(Nr,1);
            elseif strcmp(cell_flag,'mcells')
                k = ~isnan(mcells(:,1));
            elseif strcmp(cell_flag,'odorcells')
               k = mcells(:,1) > 0; 
            elseif strcmp(cell_flag,'nonodorcells')
                k = mcells(:,1) == 0;
            elseif strcmp(cell_flag,'nonmcells')
                k = isnan(mcells(:,1));
            end

            if sum(k) > 0
                % BREAK UP THE DATA CONTAINING MULTIPLE CELLS
                count = 1;
                dff = {};
                for c = 1:size(DFF{a,d},1)
                    for r = 1:size(DFF{a,d}{c},1)
                        dff{count,1} = DFF{a,d}{c}(r,:,:);
                        count = count + 1;
                    end
                end
                % ------------------
                
                dff = dff(k);
                
                for c = 1:sum(k)                                            % For each selected cell
                    dff{c} = squeeze(dff{c}(1,:,:))';                       % turn DFF to [Trials x time]
                    dff{c} = sgolayfilt(dff{c}, 1, 21, [], 2);              % Filter
    
                    [minD,maxD] = get_DFF_minmax(dff{c},time,timepoints);   % NORMALIZE ODOR DFF TO BASELINE AND COMPUTE MIN/MAX
                    minD = abs(minD);                                       % Keep absolute value of hyperpolarization
                    
                    MinD{a,d}(c,1) = nanmean(minD);                         % Get average across trials
                    MaxD{a,d}(c,1) = nanmean(maxD);                     
                    
                    FminD{a,d}(c,1) = 100*sum(~isnan(minD))/length(minD);   % Get Rate of hyperpolarizations across trials
                    
                    MD{a,d}(c,:) = mean(dff{c},1);                          % Store mean DFF over all trials
                end
            end
        end
    end
end

%% PLOT MEAN DFF ACROSS DAYS
figure;

for d = 1:maxd
    dtemp = MD(:,d); 
    dtemp = cell2mat(dtemp);
            
    subplot(maxd,4,4*(d-1)+1); hold on;
    ylabel(['Day ',num2str(d)]);
    if ~isempty(dtemp)
        for a = 1:size(dtemp,1)
            plot(time,dtemp(a,:),'Color',[0.8 0.8 0.8])
        end
        fill_plot(time,nanmean(dtemp,1),SEM(dtemp),'k')
        ylim([-1 1]); xlim([time(1) time(end)]);
    end
end

%% PLOT MINIMUM DFF ACROSS DAYS AND COMPARE NAIVE-TRAINED
A = cell(1,maxd);
for d = 1:maxd                                                          % For each day
    atemp = MinD(:,d);                                                  % Keep all mice
    atemp = cell2mat(atemp);                                            % Turn to matrix
    A{d} = atemp;                                                       % Store are cell entry
end

k = cellfun(@length,A);                                                 % Count #cells in each day
kmax = max(k);                                                          % Max # cells on a single day
Dd = nan(kmax,maxd);
for d = 1:maxd
    Dd(1:k(d),d) = A{d};
end

[Rs,Ps] = Spearman_over_days(Dd);

subplot(maxd,4,[2 3 6 7]); hold on;
for d = 1:maxd
    plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , Dd(:,d),'o','Color',[0.8 0.8 0.8]);
end
fill_plot(1:maxd,nanmean(Dd),SEM(Dd),'b')
title([num2str(Rs),' (',num2str(Ps),')']);
ylabel('|DFF MIN|');
xlim([0.5 maxd+0.5]);

subplot(maxd,4, [4 8]);
plot_naive_trained(Dd);

%% PLOT MAXIMUM DFF ACROSS DAYS AND COMPARE NAIVE-TRAINED
A = cell(1,maxd);
for d = 1:maxd                                                          % For each day
    atemp = MaxD(:,d);                                                  % Keep all mice
    atemp = cell2mat(atemp);                                            % Turn to matrix
    A{d} = atemp;                                                       % Store are cell entry
end

k = cellfun(@length,A);                                                 % Count #cells in each day
kmax = max(k);                                                          % Max # cells on a single day
Dd = nan(kmax,maxd);
for d = 1:maxd
    Dd(1:k(d),d) = A{d};
end

[Rs,Ps] = Spearman_over_days(Dd);

subplot(maxd,4,[10 11 14 15]); hold on;
for d = 1:maxd
    plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , Dd(:,d),'o','Color',[0.8 0.8 0.8]);
end
fill_plot(1:maxd,nanmean(Dd),SEM(Dd),'b')
title([num2str(Rs),' (',num2str(Ps),')']);
ylabel('DFF MAX');
xlim([0.5 maxd+0.5]);

subplot(maxd,4, [12 16]);
plot_naive_trained(Dd);

%% PLOT RATE OF HYPEROLARIZATIONS ACROSS DAYS AND COMPARE NAIVE-TRAINED (REVISIONS)
A = cell(1,maxd);
for d = 1:maxd                                                          % For each day
    atemp = FminD(:,d);                                                  % Keep all mice
    atemp = cell2mat(atemp);                                            % Turn to matrix
    A{d} = atemp;                                                       % Store are cell entry
end

k = cellfun(@length,A);                                                 % Count #cells in each day
kmax = max(k);                                                          % Max # cells on a single day
Dd = nan(kmax,maxd);
for d = 1:maxd
    Dd(1:k(d),d) = A{d};
end

[Rs,Ps] = Spearman_over_days(Dd);

subplot(maxd,4,[18 19 22 23]); hold on;
for d = 1:maxd
    plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , Dd(:,d),'o','Color',[0.8 0.8 0.8]);
end
fill_plot(1:maxd,nanmean(Dd),SEM(Dd),'b')
title([num2str(Rs),' (',num2str(Ps),')']);
ylabel('Freq. Hyperp.');
xlim([0.5 maxd+0.5]);

subplot(maxd,4, [20 24]);
plot_naive_trained(Dd);
