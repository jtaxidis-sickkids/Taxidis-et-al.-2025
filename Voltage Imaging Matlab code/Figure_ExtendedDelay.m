function Figure_ExtendedDelay(INclass)

cols = ['k';'r'];

%% FILES WITH ON-OFF TRIALS
if strcmp(INclass,'PV')
    ID{1} = 'PV1';      day{1} = '07_03_2021';  sess(1,:) = [1,2];  maxtr(1,:) = [20, 12];    
    ID{2} = 'PV1';      day{2} = '07_03_2021';  sess(2,:) = [3,4];  maxtr(2,:) = [20, 12];    
    ID{3} = 'PV1';      day{3} = '07_03_2021';  sess(3,:) = [9,10]; maxtr(3,:) = [12, 9];    
    ID{4} = 'PV1';      day{4} = '07_05_2021';  sess(4,:) = [5,6];  maxtr(4,:) = [12,8];    
    ID{5} = 'PV1';      day{5} = '07_06_2021';  sess(5,:) = [1,2];  maxtr(5,:) = [8, 8];    
    ID{6} = 'PV1';      day{6} = '07_06_2021';  sess(6,:) = [5,6];  maxtr(6,:) = [8, 8];
    
    ID{7} = 'PV2';      day{7} = '07_05_2021';  sess(7,:) = [4,5];  maxtr(7,:) = [12, 5];
    
    ID{8} = 'PV3';      day{8} = '07_05_2021';  sess(8,:) = [1,2];  maxtr(8,:) = [12, 6];
    
    ID{9} = 'PV5';      day{9} = '07_06_2021';  sess(9,:) = [1,2];  maxtr(9,:) = [12, 12];
    ID{10} = 'PV5';     day{10} = '07_06_2021'; sess(10,:) = [7,8]; maxtr(10,:) = [12, 8];    
    
elseif strcmp(INclass,'SOM')
    ID{1} = 'ASAP3-1';  day{1} = '05_18_2021';  sess(1,:) = [2,3];  maxtr(1,:) = [12, 8];    
    ID{2} = 'ASAP3-1';  day{2} = '05_18_2021';  sess(2,:) = [8,9];  maxtr(2,:) = [8, 8];
    
    ID{3} = 'ASAP3-3';  day{3} = '05_11_2021';  sess(3,:) = [6,7];  maxtr(3,:) = [12, 12];
    ID{4} = 'ASAP3-3';  day{4} = '05_14_2021';  sess(4,:) = [4,5];  maxtr(4,:) = [12, 8];
    ID{5} = 'ASAP3-3';  day{5} = '05_15_2021';  sess(5,:) = [1,2];  maxtr(5,:) = [12, 8];
    ID{6} = 'ASAP3-3';  day{6} = '05_15_2021';  sess(6,:) = [3,4];  maxtr(6,:) = [12, 8];
    ID{7} = 'ASAP3-3';  day{7} = '05_15_2021';  sess(7,:) = [6,7];  maxtr(7,:) = [8, 8];% CRAPPY
    ID{8} = 'ASAP3-3';  day{8} = '05_18_2021';  sess(8,:) = [1,3];  maxtr(8,:) = [8, 8];% Extended one crappy    
end

%% PLOT PROCESSED VIDEOS
% for d = 4%length(day)
%     VOLPY_processing(ID{d},day{d},sess(d,1))
%     VOLPY_processing(ID{d},day{d},sess(d,2))
% end


%% DETECT MCELLS IN EXTENDED DELAYS
% for d = 1:length(day)
%     Get_Modulation_2(ID{d},day{d},sess(d,2),maxtr(d,2))
% end

%% LOAD AND ORGANIZE
ld = length(day);
for del = 1:2
    for d = 1:ld
        asapfile = get_ASAPfile(ID{d},day{d},sess(d,del));
        disp(asapfile);
        [path,videoname{d,del}] = fileparts(asapfile);
        Datafile = fullfile(path,[videoname{d,del},'_data.mat']);
        modfile = fullfile(path,[videoname{d,del},'_Mcells.mat']);

        load(Datafile);
        Mc = load(modfile,'Mcells','trials','timepoints','onbins');
        
        dff = V.DFF;
        dff = dff(1,:,1:maxtr(d,del));
        DFF{d,del} = squeeze(dff)';
        
        r = V.R;
        r = r(1,:,1:maxtr(d,del));
        R{d,del} = squeeze(r)';
        
        m = B.MOT;
        M{d,del} = m(:,1:maxtr(d,del));
        
        
        trials{d,del} = Mc.trials; 
        [~ ,order{d,del}] = sort(trials{d,del});
        
        MC(d,del) = Mc.Mcells(1,1);
        F(d,del) = Mc.Mcells(1,2);
    end
    
    bins{del} = V.bins;
    time{del} = V.time;
    timepoints(del,:) = Mc.timepoints;
    onbins{del} = find(Mc.onbins);
    
    bins{del}(end) = [];    
    lt(del) = length(bins{del});
end

% R = DFF;
% bins = time;

%% NORMALIZE
for d = 1:ld
    MR = mean(R{d,1}(:,onbins{1}),1);       % Mean rate over short-delay trials over only onbins
    MR = max(MR);                        %     
    for del = 1:2
        R{d,del} = R{d,del} / MR;
    end
    
    MM = mean(M{d,1}(onbins{1},:),2);       % 
    MM = max(MM);                        %     
    for del = 1:2
        M{d,del} = M{d,del} / MM;
    end
end

%% MEAN RATES BY DELAY
Mtt = cell(ld,2);
Stt = cell(ld,2);
for d = 1:ld
    for del = 1:2
        for tt = 1:2
            r = R{d,del}(trials{d,del} == tt,:);   % Keep rate over those odors and onbins
            Mtt{d,del}(tt,:) = mean(r,1);                               
            Stt{d,del}(tt,:) = SEM(r);                               
        end
    end
end

%% GET SIGNIFICANCE BETWEEN ODOR A AND ODOR B FOR EACH DELAY
Podors = cell(1,2);
for del = 1:2
    Podors{del} = ones(ld,length(onbins{del}));
    for d = 1:ld
        for t = 1:length(onbins{del})
            R1 = R{d,del}(trials{d,del} == 1,onbins{del}(t));
            R2 = R{d,del}(trials{d,del} == 2,onbins{del}(t));           
            Podors{del}(d,t) = ranksum(squeeze(R1),squeeze(R2));
        end
        [~,Podors{del}(d,:)] = fdr(Podors{del}(d,:));
    end
end

%% PLOT ALL TRIALS EACH DAY
figure;
for d = 1:ld
    subplot(10,ld,d + [0 ld 2*ld]);    hold on;
    title([ID{d},' ',day{d},' ',videoname{d},' - ',videoname{d,2}]);

    plot_short_long_delays(R(d,:),MC(d,:),F(d,:),bins,timepoints)
end

%% PLOT FIRSTODOR-RATES SHORTvsLONG EACH DAY
for d = 1:ld
    subplot(10,ld,3*ld+d + [0 ld 2*ld]); hold on;
    
    R2 = cell(1,2);
    for del = 1:2 
        R2{del} = R{d,del}(order{d,del},:);
    end
    plot_short_long_delays(R2,MC(d,:),F(d,:),bins,timepoints)
 
    count = 0;
    for del = 1:2
        for tt = 1:2
            count = count + sum(trials{d,del} == tt);
            line([bins{1}(1) timepoints(2,6)],[1 1]*count+0.5,'color','w','linewidth',2); % Add lines between trial types
        end
    end
    
    subplot(10,ld,6*ld+[d ld+d]); hold on;
    plot_mrate_traces(Mtt{d,1},Stt{d,1},bins{1},timepoints(1,:),cols);      % COMPARE SHORT DELAY RATES
    k = (Podors{1}(d,:) < 0.05);
    plot(bins{1}(onbins{1}(k)), ones(1,sum(k))*(max(Mtt{d,1}(:))+0.2),'.k');

    subplot(10,ld,8*ld+[d ld+d]); hold on;
    plot_mrate_traces(Mtt{d,2},Stt{d,2},bins{2},timepoints(2,:),cols);      % COMPARE EXTENDED DELAY RATES
    k = (Podors{2}(d,:) < 0.05);
    plot(bins{2}(onbins{2}(k)), ones(1,sum(k))*(max(Mtt{d,2}(:))+0.2),'.k');
end


%% PLOT RATES AND FIELD SHIFTS
[~,order] = sort(F(:,1));
F = F(order,:);
MC = MC(order,:);

Mtt = cell(1,2);
for del = 1:2
    for d = 1:ld
        Mtt{del}(d,:) = mean(R{d,del},1);                               
    end
    Mtt{del} = Mtt{del}(order,:);
end

figure;
subplot(131); hold on;
plot_seq_rates(Mtt{1},[],[],bins{1},timepoints(1,:),1); % Plot trials split by type (color)
plot(F(~isnan(MC(:,1)),1),find(~isnan(MC(:,1))),'k.');                                    % Plot the fields of the long delay M-cells on the firing rates of the longer delays
plot(F(isnan(MC(:,1)),1),find(isnan(MC(:,1))),'ko');                                    % Plot the fields of the long delay M-cells on the firing rates of the longer delays
xlim([bins{1}(1) timepoints(1,5)]);

subplot(132); hold on;
plot_seq_rates(Mtt{2},[],[],bins{2},timepoints(2,:),1); % Plot trials split by type (color)
plot(F(~isnan(MC(:,2)),2),find(~isnan(MC(:,2))),'yo','MarkerSize',4,'Markerfacecolor','y');                                  
plot(F(isnan(MC(:,2)),2),find(isnan(MC(:,2))),'bo','MarkerSize',4);                                   
xlim([bins{2}(1) timepoints(2,5)]);

%% BAYESIAN DECODING OF EXTENDED DELAY BASED ON SHORT ONE

   





%%
function plot_short_long_delays(R,MC,F,bins,timepoints)
z = 0;
for del = 1:2                                                         % For each delay
    maxtr = size(R{del},1);

    imagesc(bins{del}, z + (1:maxtr), R{del},[0,1.5]);
    for x = 1:4
        line(timepoints(del,x)*[1 1], z+[0.5 maxtr+0.5], 'Color','w','linewidth',1.5); % Make lines for timepoints
    end
    line([bins{del}(1) bins{del}(end)], z*[1 1]+0.5, 'Color','w','linewidth',1.5); % Make lines for timepoints
    
    if isnan(MC(del))
        line(F(del)*[1 1] , z+[0.5 maxtr+0.5],'Color','w','Linewidth',2,'Linestyle','--');
    else
        line(F(del)*[1 1] , z+[0.5 maxtr+0.5],'Color','w','Linewidth',2,'Linestyle','-');
    end
    z = z + maxtr;
end
set(gca,'YDir','reverse');
xlim([bins{2}(1) timepoints(2,5)]);
ylim([0.5 z+0.5]);

