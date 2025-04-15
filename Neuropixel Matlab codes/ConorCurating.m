% Curating each session for Jiannis
% Read in Neuropixel files, LFP, behavior.mat files and WinEDR file
% Output one big .mat for each session that has everything

mainfolder = 'E:\MATLAB\Neuropixels Data Jiannis\Data\PreTraining';

%% Parameters

session = 'ZD028_0822';

sRate = 30000;
LFPrate = 2500;

%% Load Neuropixel data
pixelfolder = fullfile(mainfolder,session,strcat(session,'_g0_imec0'));

% cids, 2 is good, 1 is mua
[cids, cgs] = readClusterGroupsCSV(fullfile(pixelfolder,'cluster_group.tsv'));

% Times for all detected spikes
ss = readNPY(fullfile(pixelfolder,'spike_times.npy'));                % in SAMPLES
spt_all = double(ss)/sRate;                                             % converted to seconds

% Cluster ID assigned for each spike (should be same numel as spt_all)
spike_clusters = readNPY(fullfile(pixelfolder,'spike_clusters.npy'));

% Features from across all "Template matched"
% data. Cluster id can be used to access values for the SUA.
[cluster_id, cluster_ampl, cluster_chan, cluster_depth] = tsvread(fullfile(pixelfolder,'cluster_info.tsv'),...
 'cluster_id', 'Amplitude', 'ch', 'depth');

pixeldata = struct('cgs',{cgs},'cids',{cids},'cluster_ampl',{cluster_ampl},'cluster_chan',{cluster_chan},...
    'cluster_depth',{cluster_depth},'cluster_id',{cluster_id},'spt_all',{spt_all},'spike_clusters',{spike_clusters},'sRate',{sRate});

%% Load LFPs
load(fullfile(pixelfolder,strcat(session,'_g0_t0.imec0.mat')));
LFPs = dat;

[~,LFPtimepoints] = size(LFPs);
LFPtime = [1/LFPrate:1/LFPrate:LFPtimepoints/LFPrate];

%% Load WinEDR File

edrdir = struct2cell(dir(fullfile(mainfolder,session,'WinEDR files','*.edr')));                       % Directory of .edr files
[rawedr,~] = import_edr(fullfile(edrdir{2,1},edrdir{1,1}));          % Save info into rawedr cell for each ses


%% Load Behavior Files

sessiontypes = {'no light','no light','no light','no light','no light','no light','no light','no light'};
% sessiontypes = {'no light','no light','no light','no light','no light',...
%     'rebound','rebound','rebound','rebound',...
%     'onset','onset','onset','onset',...
%     'full odor','full odor','full odor'};

tempmatdir = dir(fullfile(mainfolder,session,'Behavior Matlab Files'));     % Full directory of matlab folders
matdir = struct2cell(tempmatdir([tempmatdir.isdir]));                       % Get ones that are folders

[~,numsessions] = size(matdir);
numsessions = numsessions -2;        % Because some folders '.' and '..'

allbeh = struct('licktimes',{[]},'first_lick_within_window',{[]},'Data',{[]},'progress',{[]},'sessiontype',{[]},'give_stim',{[]});

for ses=1:numsessions
    bevdir = struct2cell(dir(fullfile(matdir{2,ses+2},matdir{1,ses+2},'*.mat')));   % Directory should just have one .mat file inside
    thisbeh = load(fullfile(bevdir{2,1},bevdir{1,1}));
    allbeh.licktimes = cat(1,allbeh.licktimes,thisbeh.licktimes);
    allbeh.first_lick_within_window = cat(1,allbeh.first_lick_within_window,thisbeh.first_lick_within_window);
    allbeh.Data = cat(3,allbeh.Data,thisbeh.Data);
    allbeh.progress = cat(1,allbeh.progress,thisbeh.progress);
    
    [numtrials,~] = size(thisbeh.progress);
    thisbehsessiontype = cell(numtrials,1);
    thisbehsessiontype(:) = {sessiontypes{ses}};
    allbeh.sessiontype = cat(1,allbeh.sessiontype,thisbehsessiontype);
    
    if strcmp(sessiontypes{ses},'no light')
        allbeh.give_stim = cat(1,allbeh.give_stim,zeros(numtrials,1));
    elseif strcmp(sessiontypes{ses},'full odor')
        allbeh.give_stim = cat(1,allbeh.give_stim,ones(numtrials,1));
    elseif strcmp(sessiontypes{ses},'rebound') || strcmp(sessiontypes{ses},'onset') || strcmp(sessiontypes{ses},'hyper')
        allbeh.give_stim = cat(1,allbeh.give_stim,thisbeh.give_stim);
    end
end


%% Now detect TTL pulses and generate trialt and odort

edrpulses = [];
edrttlthresh = 2500;
for a=2:length(rawedr(:,2))-11
    % If previous timepoint is below threshold and 10 next timepoints are all above threshold
    if rawedr(a-1,2)<=edrttlthresh && sum(rawedr(a:a+9,2)>edrttlthresh)==10
        edrpulses = cat(1,edrpulses,rawedr(a,1));
    end
end

LFPpulses = [];
LFPttlthresh = 30;
for a=2:LFPtimepoints-11
    % If previous timepoint is below threshold and 10 next timepoints are all above threshold
    if LFPs(385,a-1)<=LFPttlthresh && sum(LFPs(385,a:a+9)>LFPttlthresh)==10
        LFPpulses = cat(1,LFPpulses,LFPtime(a));
    end
end

%%  Check that the pulses are the same and correct visually

% plot(LFPtime,LFPs(385,:));
% hold on;
% plot(rawedr(:,1),rawedr(:,2)/90);
% plot(rawedr(:,1),rawedr(:,4)/90);
% scatter(edrpulses,ones(369,1)*70);


% badedrpulses = [175];       % 101 for ZD032_0827 and 175 for ZD032_0828
% tempedrpulses = edrpulses;
% tempedrpulses(badedrpulses) = [];

% badLFPpulses = [3885:3896];       % Because started recordings at same time, but ended LFP late
% tempLFPpulses = LFPpulses;
% tempLFPpulses(badLFPpulses) = [];

%scatter(edrpulses,LFPpulses);
%disp(corrcoef(edrpulses,LFPpulses));


LFPpulses = tempLFPpulses;

%% Get odort

totaltrials = length(allbeh.give_stim);

odort = [];       % 1 if odor A, or 2 if odor B

odorathresh = max(rawedr(:,4))-((max(rawedr(:,4))-min(rawedr(:,4)))/4);    % Threshold is 25% of range from max value (this could be adjusted)
odorbthresh = max(rawedr(:,5))-((max(rawedr(:,5))-min(rawedr(:,5)))/4);

for a=2:length(rawedr(:,4))-2002
    % If prev timepoint is below thresh and 1999 next timepoints
    % are all above threshold and 2001th one is below threshold
    if rawedr(a-1,4)<odorathresh && sum(rawedr(a:a+1998,4)>odorathresh)==1999 && rawedr(a+2001,4)<odorathresh
        odort = cat(1,odort,[rawedr(a,1)+1,1]);
    % Now for odor B
    elseif rawedr(a-1,5)<odorbthresh && sum(rawedr(a:a+1998,5)>odorbthresh)==1999 && rawedr(a+2001,5)<odorbthresh
        odort = cat(1,odort,[rawedr(a,1)+1,2]);
    end
end

%% HERE DELETE ANY WRONG ODORS!!!!!

% badodors = [361];       % 361 for ZD032_0828
% tempodort = odort;
% tempodort(badodors,:) = [];

% scatter(odort(:,1),[1:length(odort(:,1))]);

%% Confirm that they are good now

% edrpulses = tempedrpulses;
% LFPpulses = tempLFPpulses;
% odort = tempodort;

%% Get trialt

trialt = odort(1:2:end,1) -1;       % Minus 1 because trial starts 1 sec before first odor onset


%% Now package into one big struct

fulldaydata = struct('pixeldata',{pixeldata},'allbeh',{allbeh},'LFPs',{LFPs},'LFPtime',{LFPtime},'LFPpulses',{LFPpulses},...
    'edrtime',{rawedr(:,1)},'edrpulses',{edrpulses},'rawedr',{rawedr},'trialt',{trialt},'odort',{odort});

%% Save

savefilename = strcat(session,'_curated');
save(fullfile(mainfolder,session,savefilename),'fulldaydata','-v7.3');     % If doesnt save because too big, just add '-v7.3' to end of save
