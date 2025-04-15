function Figure_ONOFF(INclass)

addpath('CircStat2012a');
addpath('DeriveLFP');

%% FILES WITH ON-OFF TRIALS
if strcmp(INclass,'PV')
    ID{1} = 'PV1';      day{1} = '07_05_2021';  sess(1) = 2;    maxtr(1) = 10;  
    ID{2} = 'PV1';      day{2} = '07_05_2021';  sess(2) = 4;    maxtr(2) = 5;    
    ID{3} = 'PV1';      day{3} = '07_05_2021';  sess(3) = 7;    maxtr(3) = 12; 
    ID{4} = 'PV1';      day{4} = '07_06_2021';  sess(4) = 3;    maxtr(4) = 12; 
    
    ID{5} = 'PV2';      day{5} = '07_05_2021';	sess(5) = 3;    maxtr(5) = 12;%8;  
    
    ID{6} = 'PV3';      day{6} = '07_05_2021';	sess(6) = 3;    maxtr(6) = 20;  
    
%     ID{7} = 'PV5';      day{7} = '07_06_2021';	sess(7) = 3;    maxtr(7) = 5;  % minfreq = 7
    ID{7} = 'PV5';      day{7} = '07_06_2021';	sess(7) = 4;    maxtr(7) = 6;  
    ID{8} = 'PV5';      day{8} = '07_06_2021';	sess(8) = 9;    maxtr(8) = 20;  
    
elseif strcmp(INclass,'SOM')
    ID{1} = 'ASAP3-3';  day{1} = '05_18_2021';  sess(1) = 2;    maxtr(1) = 10; 
    ID{2} = 'ASAP3-3';  day{2} = '05_18_2021';  sess(2) = 7;    maxtr(2) = 10;  
    ID{3} = 'ASAP3-3';  day{3} = '05_18_2021';  sess(3) = 10;   maxtr(3) = 12;  % minfreq = 1/4  
 
    ID{4} = 'ASAP3-1';  day{4} = '05_18_2021';  sess(4) = 4;    maxtr(4) = 20;  % minfreq = 1 (CRAPPY SPIKES)
    ID{5} = 'ASAP3-1';  day{5} = '05_18_2021';  sess(5) = 7;    maxtr(5) = 20;  
    ID{6} = 'ASAP3-1';  day{6} = '05_18_2021';  sess(6) = 10;   maxtr(6) = 10;  
end

ls = length(day);                                                           % Number of sessions (and cells)

%% PLOT VOLPY PROCESSING FROM EXAMPLE SESSIONS
if strcmp(INclass,'PV')
    Plot_VOLPY_processing(ID{3},day{3},sess(3));
elseif strcmp(INclass,'SOM')
    Plot_VOLPY_processing(ID{6},day{6},sess(6));
end

%% LOAD AND ORGANIZE
DFF = cell(1,ls);
R = cell(1,ls);
M = cell(1,ls);
SP = cell(1,ls);
videoname = cell(1,ls);
maxon = zeros(1,ls);

Ronoff = cell(1,ls);
Donoff = cell(1,ls);
Monoff = cell(1,ls);
MR = cell(1,ls);
SR = cell(1,ls);
MD = cell(1,ls);
SD = cell(1,ls);
MM = cell(1,ls);
SM = cell(1,ls);

for s = 1:ls
    asapfile = get_ASAPfile(ID{s},day{s},sess(s));
    disp(asapfile);
    [path,videoname{s}] = fileparts(asapfile);
    Datafile = fullfile(path,[videoname{s},'_data.mat']);
    load(Datafile,'V','B');
            
    dff = V.DFF;                                                            % Keep raw DFF
    r = V.R;                                                                % ...rates...
    sp = V.SP;                                                              % ...spike times...
    m = B.MOT;                                                              % ...and motion
    
    dff(isnan(dff) | isinf(dff)) = 0;   
    
    dff = dff(1,:,1:maxtr(s));                                              % Keep first cell (only 1 anyway) and up to max trial
    r = r(1,:,1:maxtr(s));                                          
    m = m(:,1:maxtr(s));
    sp = sp(1,1:maxtr(s));
    
    DFF{s} = squeeze(dff)';                                                 % Turn to [trials x time] and store
    R{s} = squeeze(r)';
    M{s} = m';                       
    SP{s} = sp;
end

bins = V.bins;
time = V.time;

lb = length(bins);
lt = length(time);

pars = B.pars;
timepoints = [pars.stim1_on, pars.delay_on, pars.stim2_on, pars.stim2_off, ...
    pars.lick_on, pars.vacuum_on,  pars.trial_dur];

Fs = 1000;

%% SELECT TIME TO AVERAGE OVER
onbins = (bins >= timepoints(1) & bins < timepoints(3));                    % Use bins from first odor onset to end of delay
onbins = find(onbins);

ontime = (time >= timepoints(1) & time < timepoints(3));                    % Use bins from first odor onset to end of delay
ontime = find(ontime);

%% DESPIKE THE DFF
for s = 1:ls
    Ntr = size(DFF{s},1);
    DFFconcat = DFF{s}';                                                    % [time x trials]
    DFFconcat = DFFconcat(:);                                               % Concatenate all trials in one column
    
    SPconcat = cell(1,Ntr);
    for tr = 1:Ntr
        SPconcat{tr} = SP{s}{tr} + (tr-1)*lt/Fs ;                           % Turn spiketimes from per-trial to full-time
    end
    SPconcat = cell2mat(SPconcat)';                                         % Concatenate spike times and stack in a single column

    despDFF = despike(DFFconcat,SPconcat,Fs);    
    DFF{s} = reshape(despDFF,[lt,Ntr])';                                    % Reshape to turn back to [trials x time]
end

%% SMOOTH 
for s = 1:ls
    Ntr = size(DFF{s},1);
    M{s} = sgolayfilt(M{s}, 1, 1*Fs+1, [], 2);                              % Lowpass motion with Savitzky-Golay filter
    DFF{s} = sgolayfilt(DFF{s}, 1, 21, [] , 2);                             
    for tr = 1:Ntr
        R{s}(tr,:) = smooth(R{s}(tr,:),5,'moving');
    end
end

%% NORMALIZE EACH CELL BY MAX AVERAGE TRACE
for s = 1:ls    
    DFF{s} = DFF{s} / max(mean(DFF{s}(:,ontime),2));                        % AVERAGE OVER ONBINS AND THEN TAKE MAX OVER ALL TRIALS
    R{s} = R{s} / max(mean(R{s}(:,onbins),2));           
    M{s} = M{s} / max(mean(M{s}(:,ontime),2));            
end

%% ORGANIZE PER SESSION
for s = 1:ls
    % ORGANIZE TRIAL SEQUENCE
    onoff = ones(1,maxtr(s));                                               
    onoff(2:2:maxtr(s)) = 2;                                                % Make vector with ON-OFF flags
    maxon(s) = sum(onoff == 1);                                             % Count the ON trials

    % MEAN BY ON-OFF
    MR{s} = zeros(2,lb);
    SR{s} = zeros(2,lb);
    MD{s} = zeros(2,lt);
    SD{s} = zeros(2,lt);
    MM{s} = zeros(2,lt);
    SM{s} = zeros(2,lt);
    for t = 1:2                                                             % For each type of trial (ON vs OFF)
        MR{s}(t,:) = mean(R{s}(onoff == t,:),1);                            % Keep the cell's mean firing rate over both odors
        SR{s}(t,:) = SEM(R{s}(onoff == t,:));                               % Keep the cell's SEM firing rate
        MD{s}(t,:) = mean(DFF{s}(onoff == t,:),1);                          % Same for DFF
        SD{s}(t,:) = SEM(DFF{s}(onoff == t,:));                               
        MM{s}(t,:) = mean(M{s}(onoff == t,:),1);                            % And for motion
        SM{s}(t,:) = SEM(M{s}(onoff == t,:));                              
    end

    % REORDER BY ON FOLLOWED BY OFF TRIALS
    [~,order] = sort(onoff);
    Ronoff{s} = R{s}(order,:);
    Donoff{s} = DFF{s}(order,:);
    Monoff{s} = M{s}(order,:);
end

%% GET SIGNIFICANCE
pR = ones(ls,length(onbins));
pD = ones(ls,length(ontime));
pM = ones(ls,length(ontime));    

for s = 1:ls
    onoff = ones(1,maxtr(s));                                               
    onoff(2:2:maxtr(s)) = 2;                                                % Make vector with ON-OFF flags
    
    for t = 1:length(onbins)                                                % For each modulation bin
        r1 = R{s}(onoff == 1,onbins(t));                                    % Keep rates in ON trials along that bin
        r2 = R{s}(onoff == 2,onbins(t));                                    % Same for OFF trials
        pR(s,t) = ranksum(r1,r2);                                           % Compare them
    end
    [~,pR(s,:)] = fdr(pR(s,:));                                             % FDR Correction
        
    for t = 1:length(ontime)                                                % Repeat for DFF and motion
        d1 = DFF{s}(onoff == 1,ontime(t));
        d2 = DFF{s}(onoff == 2,ontime(t));
        pD(s,t) = ranksum(d1,d2);
        
        m1 = M{s}(onoff == 1,ontime(t));
        m2 = M{s}(onoff == 2,ontime(t));
        pM(s,t) = ranksum(m1,m2);
    end
    [~,pD(s,:)] = fdr(pD(s,:));
    [~,pM(s,:)] = fdr(pM(s,:));
end

%% GET AVERAGE ON-OFF ACROSS ENTIRE TRIALS
mRon = nan(ls,1);
mRoff = nan(ls,1);
mDon = nan(ls,1);
mDoff = nan(ls,1);
mMon = nan(ls,1);
mMoff = nan(ls,1);

for s = 1:ls                                                                    
    Ron = Ronoff{s}(1:maxon(s),onbins);                                     % Keep rate over ON trials and modulation bins
    mR = mean(Ron,2);                                                       % average over bins
    mRon(s) = mean(mR,1);                                                   % and then over trials (both odors in there)
    Roff = Ronoff{s}(maxon(s)+1:end,onbins);                                % Same for OFF trials
    mR = mean(Roff,2);   
    mRoff(s) = mean(mR,1);
    
    Don = Donoff{s}(1:maxon(s),ontime);                                     % Repeat for DFF
    mD = mean(Don,2);
    mDon(s) = mean(mD,1);
    Doff = Donoff{s}(maxon(s)+1:end,ontime);
    mD = mean(Doff,2);
    mDoff(s) = mean(mD,1);
    
    Mon = Monoff{s}(1:maxon(s),ontime);                                     % Repeat for motion
    mM = mean(Mon,2);
    mMon(s) = mean(mM,1);
    Moff = Monoff{s}(maxon(s)+1:end,ontime);
    mM = mean(Moff,2);
    mMoff(s) = mean(mM,1);
end

%% PLOT ALL TRIALS EACH DAY
figure;
cols = ['k';'r'];

for s = 1:ls
    subplot(12,ls,s + [0 ls 2*ls]);
    plot_seq_rates(R{s},[],bins,timepoints,1.5);           
    title([ID{s},' ',day{s},' ',videoname{s}]);
    xlim([bins(1) timepoints(5)]);
end

%% PLOT ONvsOFF RATES
for s = 1:ls
    subplot(12,ls,3*ls+s + [0 ls 2*ls]);
    plot_seq_rates(Ronoff{s},[],bins,timepoints,1.5);                       % Plot rates in ONvsOFF trials
    line([bins(1) timepoints(7)],[1 1]*maxon(s)+0.5,'color','w','linewidth',2); % Add lines between trial types
    xlim([bins(1) timepoints(5)]);
end

%% PLOT AVERAGE ONvsOFF RATES,DFF AND MOTION PER CELL
for s = 1:ls
    subplot(12,ls,6*ls+[s ls+s]); hold on;
    plot_mrate_traces(MR{s},SR{s},bins,timepoints,cols);                    % Plot mean rates in ON vs OFF trials
    plot(bins(onbins(pR(s,:) < 0.05)), ones(sum(pR(s,:) < 0.05),1)*(max(MR{s}(:))+0.2),'.k');
    xlim([bins(1) timepoints(5)]);

    subplot(12,ls,8*ls+[s ls+s]); hold on;
    plot_mrate_traces(MD{s},SD{s},time,timepoints,cols);                    % Same for DFF
    plot(time(ontime(pD(s,:) < 0.05)), ones(sum(pD(s,:) < 0.05),1)*(max(MD{s}(:))+0.2),'.k');
    xlim([0 timepoints(5)]);

    subplot(12,ls,10*ls+s); hold on;
    plot_mrate_traces(MM{s},SM{s},time,timepoints,cols);                    % Same for motion
    plot(time(ontime(pM(s,:) < 0.05)), ones(sum(pM(s,:) < 0.05),1)*(max(MM{s}(:))+0.2),'.k');
    xlim([0 timepoints(5)]);
end

%% PLOT AVERAGE RATE IN ONvsOFF TRIALS
subplot(12,ls,11*ls+1); hold on;
for s = 1:ls
    plot(1:2,[mRon(s) mRoff(s)],'o-','Color',rand(1,3));                    % Plot mean ON vs mean OFF rates
end
[~,PR_onoff] = ttest(mRon,mRoff)                                            % Compare them

errorbar([0.9 2.1],[mean(mRon) mean(mRoff)],[SEM(mRon) SEM(mRoff)],'ok-')   % Plot means over all cells
plot_significance(PR_onoff,0.9,2.1,mean(mRon),mean(mRoff))
xlim([0.8 2.2]);
ylabel('<Rate>');

%% PLOT AVERAGE DFF IN ONvsOFF TRIALS
subplot(12,ls,11*ls+2); hold on;
for s = 1:ls
    plot(1:2,[mDon(s) mDoff(s)],'o-','Color',rand(1,3));
end
[~,PD_onoff] = ttest(mDon,mDoff)

errorbar([0.9 2.1],[mean(mDon) mean(mDoff)],[SEM(mDon) SEM(mDoff)],'ok-')
plot_significance(PD_onoff,0.9,2.1,mean(mDon),mean(mDoff))
xlim([0.8 2.2]);
ylabel('<DFF>');

%% PLOT AVERAGE MOTION IN ONvsOFF TRIALS
subplot(12,ls,11*ls+3); hold on;
for s = 1:ls
    plot(1:2,[mMon(s) mMoff(s)],'o-','Color',rand(1,3));
end
[~,PM_onoff] = ttest(mMon,mMoff)

errorbar([0.9 2.1],[mean(mMon) mean(mMoff)],[SEM(mMon) SEM(mMoff)],'ok-')
plot_significance(PM_onoff,0.9,2.1,mean(mMon),mean(mMoff))
xlim([0.8 2.2]);
ylabel('<Motion>');

%% PLOT AVERAGE DFF DIP AT ODOR ONSET
MinD_onoff = nan(ls,2);
for s = 1:ls                                                                % For each cell
    don = Donoff{s}(1:maxon(s),:); 
    [~,~,~,~,~,minD] = get_DFF_minmax(don,time,timepoints);          % NORMALIZE ODOR DFF TO BASELINE AND COMPUTE MIN/MAX
    MinD_onoff(s,1) = nanmean(minD);                                          % Get average across trials
   
    doff = Donoff{s}(maxon(s)+1:end,:); 
    [~,~,~,~,~,minD] = get_DFF_minmax(doff,time,timepoints);          % NORMALIZE ODOR DFF TO BASELINE AND COMPUTE MIN/MAX
    MinD_onoff(s,2) = nanmean(minD);                                          % Get average across trials
end

subplot(12,ls,11*ls+4); hold on;
for s = 1:ls
    plot(1:2,[MinD_onoff(s,1) MinD_onoff(s,2)],'o-','Color',rand(1,3));
end
[~,PDdip_onoff] = ttest(MinD_onoff(:,1),MinD_onoff(:,2))

errorbar([0.9 2.1],mean(MinD_onoff,1),SEM(MinD_onoff),'ok-')
plot_significance(PDdip_onoff,0.9,2.1,mean(MinD_onoff(:,1)),mean(MinD_onoff(:,2)))
xlim([0.8 2.2]);
ylabel('<DFF Dip>');

