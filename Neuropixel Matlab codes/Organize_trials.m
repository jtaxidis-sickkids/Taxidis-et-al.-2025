function Organize_trials(session,CA1ch)

Fs = 2500;
Fs_edr = 1000;

L = matfile(['Z:\Neuropixels Data Jiannis\ProcessedData\',session,'\ProcessedData.mat']);

load('All Channel Coords.mat','A')
channel_depths = A(:,2); % 

Nprobes = length(CA1ch);

%% KEEP FINAL RANGE OF CA1 CHANNELS AND UNITS
allchannels = 0:383;
chrange = nan(22,Nprobes);
for pr = 1:Nprobes
    ca1depth = channel_depths(CA1ch(pr));                                        % keep the channel's depth
    range = (channel_depths >= ca1depth - 100 & channel_depths <= ca1depth + 100); % Find channel depths within 100um from it
    
    chrange(1:sum(range),pr) = allchannels(range)';                                               % Keep the channel indexes;
end

if strcmp(session,'ZD031_0906'), chrange(22,1) = 211; end
% if strcmp(session,'CD235_0905'), chrange(22,1) = 191; end
% if strcmp(session,'ZD031_0908'), chrange(22,1) = 205; end
chrange

%% KEEP GOOD UNITS AND THEIR WAVEFORMS
CA1units = [];
allCA1spikes = {};
Nc = zeros(1,Nprobes);

S = L.pixeldata;
W = L.W;

% TURN TO CELL STRUCTURE FOR SESSIONS PROCESSED BEFORE THE 2-PROBE SYSTEM WAS ADDED IN THE CODE
if ~iscell(S.cids)                                 
    k = S.cids; S.cids = cell(1,1); S.cids{1} = k;
    k = S.cgs; S.cgs = cell(1,1); S.cgs{1} = k;
    k = S.spike_clusters; S.spike_clusters = cell(1,1); S.spike_clusters{1} = k;
    k = S.spt_all; S.spt_all = cell(1,1); S.spt_all{1} = k;
    k = S.cluster_id; S.cluster_id = cell(1,1); S.cluster_id{1} = k;
    k = S.cluster_chan; S.cluster_chan = cell(1,1); S.cluster_chan{1} = k;
    k = S.cluster_depth; S.cluster_depth = cell(1,1); S.cluster_depth{1} = k;
    k = S.cluster_ampl; S.cluster_ampl = cell(1,1); S.cluster_ampl{1} = k;
    k = W; W = cell(1,1); W{1} = k;
end

for pr = 1:Nprobes
    goodunits = S.cids{pr}(S.cgs{pr} == 2);                                             % Find cluster IDs with good units
    
    N = length(goodunits);
    spikes = cell(N,1);                                                         % Make cell array for spikes of units stacked by cgs
    channels = zeros(N,3);                                                      % Make matrix for their channel number/depth also stacked by cgs
    channels(:,1) = goodunits';                                                 % Store good unit IDs as column
    
    for i = 1:N                                                                 % For each good unit
        k = (S.spike_clusters{pr} == goodunits(i));                                 % Find where it appears in spike_clusters
        spikes{i} = S.spt_all{pr}(k);                                               % Keep its corresponding spike times
        
        k = find(S.cluster_id{pr} == goodunits(i));                                 % Find where the unit appears in the list of all units
        channels(i,2) = S.cluster_chan{pr}(k);                                      % Keep the corresponding channel
        channels(i,3) = S.cluster_depth{pr}(k);                                     % And its depth
        channels(i,4) = S.cluster_ampl{pr}(k);                                      % And spike amplitude
    end
    
    ca1_chan = (channels(:,2) >= chrange(1,pr) & channels(:,2) <= chrange(end,pr));          % Find units in channels within the range
    CA1units = [CA1units; channels(ca1_chan,:)];                                                   % Keep corresponding IDs/channels/depths etc
    allCA1spikes = [allCA1spikes; spikes(ca1_chan)];                                                   % And their corresponding spikes
    
    Nc(pr) = sum(ca1_chan);
    if Nc == 0
        error('No units in CA1 pyr layer!');
    end

    W{pr} = W{pr}(S.cgs{pr} == 2,:,:);                                                      % Keep the good unit waveforms (stacked in same arrangement as cgs)
    W{pr} = W{pr}(ca1_chan,:,:);                                                               % And their corresponding waveforms
end
Nc

W = cell2mat(W');

clear S

%% LOAD LFPs
LFPs = L.LFPs;

lt = size(LFPs,2);
LFPtime = 0 : 1/Fs : (lt-1)/Fs;

LFPttl = LFPs(385,:,1);
LFPs = LFPs(chrange(:,pr)+1,:,:);                                     % Keep corresponding LFP channels (+1 since channels start from 0)   

%% ORGANIZE BEHAVIOR AND GET OPTO PROTOCOL PER TRIAL
B = L.allbeh;
progress = B.progress(:,1:2);

Ntr = size(progress,1);

opto = nan(Ntr,1);
k = contains(B.sessiontype,'no light');
opto(k) = 0;
k = contains(B.sessiontype,'onset');
opto(k) = 1;
k = contains(B.sessiontype,'hyper');
opto(k) = 2;
k = contains(B.sessiontype,'rebound');
opto(k) = 3;
k = contains(B.sessiontype,'full odor');
opto(k) = 4;

opto = [opto, B.give_stim];

%% GET WINEDR TRIAL ONSET POINTS 
rawedr = L.rawedr;
edrtime = rawedr(:,1);

bothodors = rawedr(:,4)+rawedr(:,5);                            % Add the two odor signal on WinEDR

if strcmp(session,'CD234_0825')
    bothodors(4.6e6 : 4.8e6) = 150;
    bothodors(7.63e6 : 7.72e6) = 150;
end
if strcmp(session,'ZD032_0828')
    bothodors(4.03e6 : 4.04e6) = 140;
end
if strcmp(session,'ZD028_0904')
    bothodors(1 : 2.4e5) = 140;
end

cross = (bothodors > 1000);                                          % Find when the odors are ON
edr_on = find((cross - circshift(cross,1)) == 1);                    % Find when they turn ON
edr_on = edr_on(1:2:end);                                         % Skip every second odor to keep trial initiations only
edr_t_on = edrtime(edr_on);                                         % Get WinEDR timepoint of trial onset
        
eNtr = length(edr_on);                                               % Number of EDR trials
[eNtr Ntr]

if eNtr ~= Ntr
    error('EDR trials not detected propery');
end   

LFP_on = nan(Ntr,1);
LFP_t_on = nan(Ntr,1);

%% GET LFP AND WINEDR TTL PULSES ONSET
edrttl = rawedr(:,2);                                                       % Keep the WinEDR TTL signal

if strcmp(session,'CD234_0825')
    edrttl(4.6e6 : 7.715e6) = 0;
end
if strcmp(session,'ZD032_0827')
    edrttl(2.234e6 : 2.236e6) = 0;
    LFPttl(5.5878e6 : 5.5882e6) = 0;
end
if strcmp(session,'ZD032_0828')
    edrttl(4.03e6 : 4.04e6) = 0;
    LFPttl(1.009e7 : 1.01e7) = 0;
end
if strcmp(session,'CD234_0825')
    LFPttl(4.56e5 : 4.7e5) = 0;
    LFPttl(5.5e6 : 5.51e6) = 0;
    LFPttl(1.934e6 : 1.944e6) = 0;
    LFPttl(1.15e7 : 1.93e7) = 0;
end
if strcmp(session,'ZD028_0904')
    edrttl(1 : 2e5) = 0;
    LFPttl(1 : 6e5) = 0;
end
if strcmp(session,'ZD031_0906')
    edrttl(5.42e5 : 5.56e5) = 0;
    LFPttl(1.36e6: 1.37e6) = 0;
end
if strcmp(session,'ZD031_0908')
    edrttl(5.2e5 : 5.6e5) = 0;
    LFPttl(1.3e6: 1.4e6) = 0;
end
cross = (edrttl > 4000);                                                    % Find TTLs
edrttl_on = find((cross - circshift(cross,1)) == 1);                        % Find when they turn ON

cross = (LFPttl' > 10);                                                     % Find TTLs in Neuropixels
lfpttl_on = find((cross - circshift(cross,1)) == 1);                        % Find when they turn ON

if strcmp(session,'ZD028_0822')                                             % In problematic session 
    dttl = length(lfpttl_on)-length(edrttl_on);
    lfpttl_on(1 : dttl) = [];                                               % Remove extra Neuropixel TTLs from the beginning
end

lfp_on = lfpttl_on - 1*Fs;                                                  % Find point of trial onset (odor - 1 sec; works only for no-opto)

[length(edrttl_on) length(lfpttl_on)]

%% GET LFP TRIAL ONSET TIMES FOR NO-OPTO RECORDINGS
if all(opto(:,1) == 0)                                          % If no opto was given
    if strcmp(session,'ZD028_0822')
        DT = LFPtime(lfpttl_on)' - edrtime(edrttl_on);
        for tr = 1:Ntr
            [~,closest_pulse] = min(abs(edrttl_on - edr_on(tr)));            % Find the closest WinEDR TTL to trial onset point 
            lfptime_shifted = LFPtime - DT(closest_pulse);             % Shift LFP time to match the corresponding WinEDR TTL times
            [~,k] = min(abs(lfptime_shifted - edr_t_on(tr))); % Find the closest entry to the the corresponding WinEDR trial onset time
            
            k = k - 0.5*Fs;     % I dont know why but this is necessary for alignment to be (roughly) correct
            
            LFP_on(tr) = k;
            LFP_t_on(tr) = LFPtime(k);
        end
    else
        LFP_on = lfp_on;
        LFP_t_on = LFPtime(LFP_on)';
    end
end

%% GET LFP TRIAL ONSET TIMES FOR OPTO RECORDINGS (NO-LIGHT TRIALS MISSING LFP TTL PULSE)
if any(opto(:,1) > 0)    
    Dt = LFPtime(lfpttl_on)' - edrtime(edrttl_on); % Get time distances between the two TTL signals (missing lightOFF trials and with 2 pulses for lightON trials)
    k1 = find(opto(:,1) > 0,1,'first');  % find first trial with opto
    k2 = find(opto(:,1) > 0,1,'last');  % find last trial with opto
    Dt(k1+1:2:end-(Ntr-k2)) = []; % delete the second pulse (second odor) 
                                  % ONLY WORKS IF NOLIGHT IS GIVEN ONLY IN THE BEGINNIG TRIALS 
                                  % AND IN THE FINAL TRIALS (BUT NOT IN BETWEEN!)   
    DT = nan(Ntr,1);
    k = (opto(:,2) == 1 | opto(:,1) == 0);    % Find trials with lightON or without opto protocol
    DT(k) = Dt;                                 % Keep their time distance in the corresponding entries
    
    for i = find(~k)                            % For each lightOFF trial
        z = find(~isnan(DT(1:i)),1,'last');     % Find the closest preceding lightON trial
        DT(i) = DT(z);                          % Use that time distance
    end
    
    for tr = 1:Ntr                                     % For each trial
        lfptime_shifted = LFPtime - DT(tr);             % Shift LFP time to match the corresponding WinEDR TTL times
        [~,k] = min(abs(lfptime_shifted - edr_t_on(tr))); % Find the closest entry to the the corresponding WinEDR trial onset time
        LFP_on(tr) = k;
        LFP_t_on(tr) = LFPtime(k);
    end
end

%% PLOT TRIAL ONSET FROM WINEDR AND LFPS
figure;
plot(edrtime,bothodors/100)
hold on
plot(LFPtime,LFPttl)
plot(LFP_t_on,ones(Ntr,1)*20,'vk')
plot(edr_t_on,ones(Ntr,1)*21,'vb')
ylim([0 30]);

%% SPLIT SPIKES INTO TRIALS
CA1spikes = cell(sum(Nc),Ntr);
for c = 1:sum(Nc)                                                                % For each unit
    csp = allCA1spikes{c};                                                  % Keep its spikes
    for tr = 1:Ntr                                                          % For each trial
        k = (csp >= LFP_t_on(tr) & csp <= LFP_t_on(tr) + 11);               % find spikes within 11 sec from trial onset time
        CA1spikes{c,tr} = csp(k) - LFP_t_on(tr);                            % Store them with onset time being = 0
    end
end

%% COMPUTE RATES
binlen = 0.01;                                                              % 10msec rate bins
bins = 0:binlen:11;

R = nan(sum(Nc),length(bins)-1,Ntr);

for tr = 1:Ntr
    for c = 1:sum(Nc)
        R(c,:,tr) = histcounts(CA1spikes{c,tr},bins);
    end
end
R = R/binlen;

%% SPLIT LFP INTO TRIALS
ltr = 11*Fs;

L = nan(22,ltr,Ntr,Nprobes);
for tr = 1:Ntr
    L(:,:,tr,:) = LFPs(:, LFP_on(tr) + (0:ltr-1), :);
end

for pr = 1:Nprobes
    for chan = 1:10:22
        figure; hold on;
        for tr = 1:Ntr
            plot(0:(1/Fs):(11-1/Fs),L(chan,:,tr,pr)/max(L(chan,:,tr,pr)) + tr,'Color',[.5 .5 .5]);
        end
        xlim([0 3])
        ylim([0 Ntr+1])
    end
end

L = int16(L);

%% SPLIT LOCOMOTION INTO TRIALS
ltr = 11*Fs_edr;
v = rawedr(:,3);
V = nan(Ntr,ltr);
for tr = 1:Ntr
    V(tr,:) = v(edr_on(tr) + [0:ltr-1]);
end

%% SAVE
dirname = ['../ProcessedData/',session];
if ~exist(dirname,'dir')                                           
    mkdir(dirname);
end
save(fullfile(dirname,'CA1data.mat'),'V','L','R','W','CA1spikes','CA1units','Nc','progress','opto','CA1ch','chrange','bins')
