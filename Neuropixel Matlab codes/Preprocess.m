function Preprocess(session,sessiontypes,Nprobes)

%% Folders
dirname = ['Z:\Neuropixels Data Jiannis\ProcessedData\',session];
if ~exist(dirname,'dir')                                           
    mkdir(dirname);
end
datafile = fullfile(dirname,'ProcessedData.mat');

sessfolder = fullfile('Z:\Neuropixels Data Jiannis\Data',session);
newropixfolder = cell(1,Nprobes);
for pr = 1:Nprobes 
    newropixfolder{pr} = fullfile(sessfolder,[session,'_g0_imec',num2str(pr-1)]);
end

%% Load LFPs and save
LFPs = cell(1,1,Nprobes);
for pr = 1:Nprobes 
    L = matfile(fullfile(newropixfolder{pr},[session,'_g0_t0.imec',num2str(pr-1),'.mat']));
    LFPs{pr} = L.dat;
end

lt = cellfun(@(x) size(x,2),LFPs)
if Nprobes == 2 && lt(1) ~= lt(2)
   [minlt,k] = min(lt);
   LFPs{setdiff(1:2,k)}(:,minlt+1:end) = [];
end

LFPs = cell2mat(LFPs);

save(datafile,'LFPs','-v7.3'); 
clear LFPs 

%% Load Neuropixel data and waveforms and save
cids = cell(1,Nprobes);
cgs = cell(1,Nprobes);
spt_all = cell(1,Nprobes);
spike_clusters = cell(1,Nprobes);
cluster_id = cell(1,Nprobes);
cluster_ampl = cell(1,Nprobes);
cluster_chan = cell(1,Nprobes);
cluster_depth = cell(1,Nprobes);
W = cell(1,Nprobes);

for pr = 1:Nprobes
    % cids, 2 is good, 1 is mua
    [cids{pr}, cgs{pr}] = readClusterGroupsCSV(fullfile(newropixfolder{pr},'cluster_group.tsv'));
    
    % Times for all detected spikes
    spt_all{pr} = readNPY(fullfile(newropixfolder{pr},'spike_times.npy'));                % in SAMPLES
    spt_all{pr} = double(spt_all{pr})/30e3;                                             % converted to seconds
    
    % Cluster ID assigned for each spike (should be same numel as spt_all)
    spike_clusters{pr} = readNPY(fullfile(newropixfolder{pr},'spike_clusters.npy'));
    
    % Features from across all "Template matched" data.
    [cluster_id{pr}, cluster_ampl{pr}, cluster_chan{pr}, cluster_depth{pr}] = tsvread(fullfile(newropixfolder{pr},'cluster_info.tsv'),...
        'cluster_id', 'Amplitude', 'ch', 'depth');
  
    % Load waveforms of all cells
    W{pr} = nan(length(cids{pr}),40,200);                                        % Use it to get number of trials to set size
    for c = 1:length(cids{pr})                                                      % For each unit
        ww = h5info(fullfile(newropixfolder{pr}, ['Waveforms_',session,'_Imec',num2str(pr-1),'.h5']),sprintf('/%d', c-1));  % Load waveform numbered by unit index in cids (starts from 0)
    
        if size(ww,3) < 200
            ww = cat(3,ww,nan(size(ww,1),40, 200 - size(ww,3)));
        end
        
        cid = cids{pr}(c);                                                          % Keep the unit's ID number
        k = find(cluster_id{pr} == cid);                                            % Find where it is on cluster list
        maxchannel = cluster_chan{pr}(k);                                           % Find the corresponding channel
        W{pr}(c,:,:) = ww(maxchannel+1, :, :);                                      % Keep that channel (+1 cause channels start at 0) x trials x timepoints
        %     figure;hold on;for i =1:200,plot(W{pr}(c,:,i),'Color',[.8 .8 .8]);end; plot(mean(W{pr}(c,:,:),3),'k','LineWidth',2);title(num2str(maxchannel));
    end
end
pixeldata = struct('cgs',{cgs},'cids',{cids},'cluster_ampl',{cluster_ampl},'cluster_chan',{cluster_chan},...
        'cluster_depth',{cluster_depth},'cluster_id',{cluster_id},'spt_all',{spt_all},'spike_clusters',{spike_clusters});

save(datafile,'pixeldata','W','-append'); 
clear W pixeldata

%% Load WinEDR File
edrfile = dir(fullfile(sessfolder,'*.edr'));                      % Directory of .edr files
rawedr = import_edr(fullfile(sessfolder,edrfile.name));          % Save info into rawedr cell for each ses

save(datafile,'rawedr','-append'); 
clear rawedr

%% Load Behavior Files
matfiles = dir(fullfile(sessfolder,'*.mat'));     % Full directory of matlab folders
matdates = [matfiles.datenum]';
[~,order] = sort(matdates);
matfiles = {matfiles.name}';                      % Get ones that are folders
matfiles = matfiles(order)

lm = length(matfiles);

allbeh = struct('licktimes',{[]},'first_lick_within_window',{[]},'Data',{[]},'progress',{[]},'sessiontype',{[]},'give_stim',{[]});

for m = 1:lm
    mtemp = load(fullfile(sessfolder,matfiles{m}));
    allbeh.licktimes = cat(1,allbeh.licktimes,mtemp.licktimes);
    allbeh.first_lick_within_window = cat(1,allbeh.first_lick_within_window,mtemp.first_lick_within_window);
    allbeh.Data = cat(3,allbeh.Data,mtemp.Data);
    allbeh.progress = cat(1,allbeh.progress,mtemp.progress);
    
    Ntr = size(mtemp.progress,1);
    opto = cell(Ntr,1);
    opto(:) = {sessiontypes{m}};
    allbeh.sessiontype = cat(1,allbeh.sessiontype,opto);

    if strcmp(sessiontypes{m},'no light')
        allbeh.give_stim = cat(1,allbeh.give_stim,zeros(Ntr,1));
    elseif strcmp(sessiontypes{m},'full odor')
        allbeh.give_stim = cat(1,allbeh.give_stim,ones(Ntr,1));
    elseif strcmp(sessiontypes{m},'rebound') || strcmp(sessiontypes{m},'onset') || strcmp(sessiontypes{m},'hyper')
        allbeh.give_stim = cat(1,allbeh.give_stim,mtemp.give_stim);
    end
end

save(datafile,'allbeh','-append'); 
