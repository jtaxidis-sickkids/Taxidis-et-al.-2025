function [MRonoff,bins] = Figure_fine_timescale_ONOFF(INclass)

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

%% LOAD AND ORGANIZE
SP = cell(1,ls);

for s = 1:ls
    asapfile = get_ASAPfile(ID{s},day{s},sess(s));
    disp(asapfile);
    [path,videoname] = fileparts(asapfile);
    Datafile = fullfile(path,[videoname,'_data.mat']);
    load(Datafile,'V','B');
    
    sp = V.SP;                                                              % ...spike times..
    sp = sp(1,1:maxtr(s));
    SP{s} = sp;
end

bins = V.bins;
time = V.time;

lb = length(bins);
lt = length(time);

pars = B.pars;
timepoints = [pars.stim1_on, pars.delay_on, pars.stim2_on, pars.stim2_off, ...
    pars.lick_on, pars.vacuum_on,  pars.trial_dur];

%% SELECT TIME TO AVERAGE OVER
binlen = 0.005;                                                             % SHORT BIN LENGTH OF 5ms
bins = 0 : binlen : 10;%timepoints(6);
lb = length(bins);

R = cell(ls,1);
for s = 1:ls                                                                % For each cell
    Ntr = length(SP{s});
    R{s} = nan(Ntr,lb-1);
    
    for tr = 1:Ntr
        h = histcounts(SP{s}{tr},bins) / binlen;                            % COMPUTE FIRING RATE
        R{s}(tr,:) = smooth(h,3);                                           % AND SMOOTH
    end
end

%% NORMALIZE
% onbins = (bins >= timepoints(1) & bins < timepoints(3));                    % Use bins from first odor onset to end of delay
% onbins = find(onbins);
% for s = 1:ls
%     R{s} = R{s} / max(mean(R{s}(:,onbins),2));
% end

%% SPLIT TO ON OFF
Ronoff = cell(ls,2);
for s = 1:ls
    % ORGANIZE TRIAL SEQUENCE
    onoff = ones(1,maxtr(s));
    onoff(2:2:maxtr(s)) = 2;                                                % Make vector with ON-OFF flags
    
    % SPLIT BY ON-OFF
    for t = 1:2                                                             % For each type of trial (ON vs OFF)
        Ronoff{s,t} = R{s}(onoff == t,:);                                   % Keep the cell's mean firing rate over both odors
    end
end

%% GET MEAN SPIKE HISTOGRAMS
MRonoff = cell(1,2);
for tt = 1:2                                                           % For each trial type
    h = cellfun(@(x) mean(x,1), Ronoff(:,tt),'UniformOutput',false);      % GET CELL'S MEAN RATE OVER ALL TRIALS OF THAT TYPE
    MRonoff{tt} = cell2mat(h);
end

%% PLOT AVERAGE FIRING RATES OF ODOR-SPECIFIC AND NON-SPECIFIC MCELLS
figure;hold on;
cols = ['k';'r'];
for tt = 1:2                                              % For each trial type (again)
    histogram('BinEdges',bins,'BinCounts',mean(MRonoff{tt},1),...
        'EdgeColor','none','FaceColor',cols(tt,:));
    fill_plot(bins(1:end-1),mean(MRonoff{tt},1),SEM(MRonoff{tt}),cols(tt,:)); % Plot mean+SEM curves
  
    axis tight;
    xlim([0.8 2.2]);
end
