
%% DAYS WITH VALVE IN HAND
%       day        video maxtrial valveinhandtrials
days = {'03_30_2022', 4,    11,     [1:5];         % no dip
    '03_30_2022', 5,    10,     [1:5];         % no dip
    '03_30_2022', 7,    9,      [1:5];         % no dip
    '03_30_2022', 9,    8,      [5:8];         % no dip
    '04_01_2022', 4,    12,     [9:12];
    '04_13_2022', 2,    11,     [5:9];
    '04_13_2022', 3,    12,     [10:12];
    '04_13_2022', 8,    12,     [5:8]};

ld = size(days,1);

%% ASAP3-4 EXAMPLES
for d = [1 8]% 1:ld
    Plot_VOLPY_processing('ASAP3-4',days{d,1},days{d,2});         % no dip
end

%% COMPARE VALVE IN HAND VS NOT
Fs = 1000;
time = 0 : 1/Fs : (11*Fs-1)/Fs;
DFF = cell(ld,2);

for d = 1:ld
    [asapfile,~] = get_ASAPfile('ASAP3-4',days{d,1},days{d,2});
    [path,videoname] = fileparts(asapfile);
    
    datafile = fullfile(path,[videoname,'_data.mat']);
    load(datafile);
    modfile = fullfile(path,[videoname,'_Mcells.mat']);
    load(modfile,'timepoints');
    
    dff = squeeze(V.DFF(1,:,:))';                                           % Keep DFF NOT DESPIKED [Trials x time]
    dff = dff(1:days{d,3},:);                                               % Keep up to max trial
    dff = sgolayfilt(dff, 1, 21, [], 2);                                    % Filter
    
    [Ntr,lt] = size(dff);
    
    DFF{d,1} = dff(setdiff(1:Ntr,days{d,4}),:);
    DFF{d,2} = dff(days{d,4},:);
end

%% PLOT TRIALS 
figure;
cols = cell_colors;
for d = 1:ld
    % COMPUTE MEAN BY FIRSTODOR
    MD = zeros(2,lt);
    SD = zeros(2,lt);
    for tt = 1:2                                                             % For each first odor
        MD(tt,:) = mean(DFF{d,tt},1);                              % Keep the cell's mean firing rate
        SD(tt,:) = SEM(DFF{d,tt});                                 % Keep the cell's SEM
    end
    SD(isnan(SD) | isinf(SD)) = 0;
    
    % PLOT EACH TRIAL TYPE
    subplot(ld,2,2*d-1); hold on;
    for tt = 1:2
        for tr = 1:size(DFF{d,tt},1)
            plot(time,DFF{d,tt}(tr,:) + 0.8*(tt-1),'Color',cols(tt,:)); % Plot trials split by type (color)
        end
    end
    
    subplot(ld,2,2*d); hold on;
    plot_mrate_traces(MD,SD,time,timepoints,cols);
end

%% COMPARE DIPS
MinD = nan(ld,2);
for d = 1:ld                                                               % For each cell
    for tt = 1:2
        [~,~,~,~,~,minD] = get_DFF_minmax(DFF{d,tt},time,timepoints);                   % NORMALIZE ODOR DFF TO BASELINE AND COMPUTE MIN/MAX
        MinD(d,tt) = nanmean(minD);                                            % Get average dip across trials
    end
end

[~,pD] = ttest(MinD(:,1),MinD(:,2))

figure;
hold on
for d = 1:ld
    plot(1:2,MinD(d,:),'o-','Color',[0.8 0.8 0.8]);                    % Plot mean ON vs mean OFF rates
end
errorbar([0.9 2.1],nanmean(MinD,1),SEM(MinD),'ok-')   % Plot means over all cells
plot_significance(pD,0.9,2.1,nanmean(MinD(:,1)),nanmean(MinD(:,2)))
xlim([0.8 2.2]);
ylabel('Average Min');
set(gca,'XTick',[1 2],'XTickLabel',{'Valve on','Valve in hand'});

