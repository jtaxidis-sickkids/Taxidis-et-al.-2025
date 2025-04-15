function Unit_analysis(ID,day,session)

addpath('CircStat2012a');
addpath('DeriveLFP');

%% PARAMETERS
freq = [0.5 3;        % Delta (based on power spectra early peak)
        4   10;       % Theta frequency range (Bittner et al 2015 but also the power spectra peaks).
        15  30;      % Beta (based on literature and mean ISI distributions)
        40  90];     % Gamma

lf = size(freq,1);
Fs = 1000;
for f = 1:lf
    [filtb{f},filta{f}] = butter(3,freq(f,:)/(Fs/2),'bandpass');                        % construct the bandpass filter
end

window = 60;

wpostmax = 50;
burstDT = 0.010;  % 10 ms

%% LOAD FILES
asapfile = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_data.mat']);
load(Datafile,'V','B');

bins = V.bins;
time = V.time;
pars = B.pars;

lb = length(bins) - 1;

%% POOL TRIALS AND KEEP SELECTED ROI AND UP TO MAX TRIAL
[DFF,SP,R,MOT,trials,outcomes,progress,licks,roi] = load_and_pool(ID,day,session);                               % Pool trials of the same ROI

DFF(isinf(DFF)|isnan(DFF)) = 0;                                             % Remove Nans and Infs
[Nr,lt,Ntr] = size(DFF);

Template = V.Template;
ROI = V.ROI(roi,:,:);
CC = V.CC(roi);

clear V B

%% CONCATENATE TRIALS AND SPIKES FOR SPIKE SUBTRACTION AND WAVEFORM ESTIMATION 
DFFconcat = reshape(permute(DFF,[2,3,1]),lt*Ntr,[]);                        % Concatenate all trials in one column

SPconcat = cell(size(SP));
for r = 1:Nr
    for tr = 1:Ntr
        SPconcat{r,tr} = SP{r,tr} + (tr-1)*lt/Fs ;                          % Turn spiketimes from per-trial to full-time 
    end
    SPconcat{r,1} = cell2mat(SPconcat(r,:));                                % Concatenate spike times and keep them in the first cell column
end
SPconcat = SPconcat(:,1);                                                   % Keep the first cell column with all spikes

%% COMPUTE ISI
ISI = cell(Nr,Ntr);

for r = 1:Nr
    for tr = 1:Ntr
        ISI{r,tr} = diff(SP{r,tr});                                         
    end
end

%% COMPUTE SPIKE CHARACTERISTICS
PeriSp = cell(Nr,1);
Trough = nan(Nr,1);
PTdist = nan(Nr,1);
Halfwidth = nan(Nr,1);
Points = nan(Nr,3);
Nsp = nan(Nr,1);

% CONTATENATE SPIKES AND DFF
for r = 1:Nr
    spconcat = round(SPconcat{r}*Fs+1);
    
    Nsp(r) = length(spconcat);
    perispike = nan(Nsp(r),2*window+1);
    
    for sp = 1:Nsp(r)
        spikewindow = spconcat(sp) + (-window:window);
        spikewindow(spikewindow < 1) = 1;
        spikewindow(spikewindow > lt*Ntr) = lt*Ntr;
        perispike(sp,:) = DFFconcat(spikewindow,r);
    end
    mperi = mean(perispike,1);
    
    mperism = smooth(mperi,5,'moving');
    mperism = -(mperism);
    [~,k] = findpeaks(mperism(window+2 + (0:wpostmax)),'SortStr','descend');   % Get peaks in descending order and their times
    if isempty(k), k = wpostmax; end
    k = k(1) + window + 2;
    trough = mperi(k) - mean(mperi(1 : window-10));    % COMPUTE TROUGH RELATIVE TO AVERAGE PRESPIKE BASELINE                                          % Waveform trough (third peak)
    ptdist = (k-window-2)/Fs;                         % Peak - Trough distance
    
    peak =  mperi(window+2);
    [~,z] = findpeaks(-abs(mperi - peak/2),'SortStr','descend');     % Find times closest to the half-peak points (essentially in their original order)
    z = sort(z(1:2));                                                    % Keep the two closest
    halfwidth = abs(z(2) - z(1))/Fs ;                 % Halfwidth in ms
    
    PeriSp{r} = perispike;
    Trough(r) = trough;
    PTdist(r) = ptdist;
    Halfwidth(r) = halfwidth;
    Points(r,:) = [k z];
end

%% COMPUTE MEAN FIRING RATE, BURST INDEX, CSI
mR = nan(Nr,1);
BI = nan(Nr,1);
CSI = nan(Nr,1);

for r = 1:Nr
    mr = squeeze(R(r,:,:));
    mR(r) = mean(mr(:));
    
    isi = ISI(r,:);
    isi = cell2mat(isi);
    bursts = [0, isi < burstDT];                                            % Get bursty spikes
    BI(r) = sum(bursts) / Nsp(r) * 100;                                     % Burst index = percentage of spikes closer than 10 ms
    
    amps = PeriSp{r}(:,window+1);                                           % Get the spike amplitudes of all spikes
    mr = find(bursts - circshift(bursts,-1) == -1);                         % Find zeros followed by a 1 (1st spike in burst)
    bf = find(bursts - circshift(bursts,-1) == 1);                          % Find ones followed by a 0 (last spike in burst)
    decr = zeros(1,length(mr));
    incr = zeros(1,length(mr));
    for b = 1:length(mr)                                                    % For each burst
        bamps = amps(mr(b) : bf(b));                                        % Keep spike amplitudes in the burst
        dbamps = diff(bamps);                                               % Compute change in amplitude relative to predecessor (excludes 1st spike in burst)
        decr(b) = sum(dbamps < 0);                                          % Count amplitude drops
        incr(b) = sum(dbamps > 0);                                          % And increases
    end
    CSI(r) = (sum(decr) - sum(incr)) / sum(bf-mr) * 100;                    % Ratio of total number of amplitude drops - increases over total number of bursty spikes
end

%% DESPIKE THE DFF FOR SPECTRAL ANALYSIS
despDFF = nan(size(DFFconcat));
for r = 1:Nr
    [despDFF(:,r),params] = despike(DFFconcat(:,r),SPconcat{r},Fs);
end

despDFF = despDFF';                                                         % Turn back to [Nr x time]

%% GET DFF OSCILLATION AMPLITUDES AND PHASES
Ampl = nan([size(despDFF),lf]);
Phase = nan([size(despDFF),lf]);

for f = 1:lf                                                                % For each frequency
    for r = 1:Nr                                                            % For each ROI
        DFFf = filtfilt(filtb{f},filta{f},despDFF(r,:));                    % bandpass all trials
        H = hilbert(DFFf);                                                  % Get hilbert transform
        
        Ampl(r,:,f) = abs(H);                                               % Get amplitudes
        Phase(r,:,f) = angle(H);                                            % Get phases
        
        Ampl(r,:,f) = sgolayfilt(Ampl(r,:,f), 1, 1*Fs+1);                   % Lowpass  amplitude signal with Savitzky-Golay filter
    end
end

%% RESHAPE DESPIKED DFF, AMPLITUDES AND PHASES INTO TRIALS
despDFF = reshape(despDFF,size(DFF));                                       % Turm to [Nr,lt,Ntr]
Ampl = reshape(Ampl,[size(DFF),lf]);                                        % Turm to [Nr,lt,Ntr,lf]
Phase = reshape(Phase,[size(DFF),lf]);

%% GET PHASES OF SPIKES
SP_ph = cell(Nr,Ntr);
pvalR = nan(Nr,lf);
pvalV = nan(Nr,lf);
mDir = nan(Nr,lf);
lVec = nan(Nr,lf);
Skew = nan(Nr,lf);

for r = 1:Nr                                                                % For each roi
    for tr = 1:Ntr
        sp = SP{r,tr}*Fs + 1;                                               % Turn its spiketimes to time-axis indexes
        sp = sp + 1;                                                        % VOLPY SPIKES ARE DETECTED AT SPIKE ONSET - 1 TIMEPOINT BEFORE PEAK
        sp(sp > lt) = lt;                                                   % make sure all spiketimes are within limits
        sp = round(sp);                                                     % make sure they are round numbers
        
        spph = squeeze(Phase(r,sp,tr,:));                                   % Keep the corresponding phases [spikes x frequencies]
        SP_ph{r,tr} = spph;                                                 % Store it
    end
    
    sp = SP_ph(r,:);
    sp = sp';
    sp = cell2mat(sp);                                                      % Concatenated spikes x frequencies
    
    for f = 1:lf
        mDir(r,f) = circ_mean(sp(:,f));                                     % Mean theta phase
        lVec(r,f) = circ_r(sp(:,f));                                        % Length of mean vector
        pvalR(r,f) = circ_rtest(sp(:,f));                                   % Rayleigh test for non-uniformity of circular data.
        pvalV(r,f) = circ_vtest(sp(:,f), 0);                                % V-test (KUIPER'S TEST) of non-uniformity with 0-degrees mean direction
        Skew(r,f) = circ_skewness(sp(:,f));                                 % circular skewness
    end
end

%% COMPUTE SPECTROGRAM
f = 0:0.1:100;
PSD = nan(Nr,4001,103,Ntr);
for r = 1:Nr
    for tr = 1:Ntr
        [~,~,t,psd] = spectrogram(despDFF(r,:,tr), 512, 256, f, Fs,'yaxis');
        
        % Interpolate
        tq = t(1):0.1:t(end);
        fq = f(1):0.025:f(end);
        [T,Fint] = meshgrid(t,f);
        [Tq,Fq] = meshgrid(tq,fq);
        psdq = interp2(T,Fint,psd,Tq,Fq);
        
        % Smooth
        psdq = imgaussfilt(psdq,'FilterSize',5);
        PSD(r,:,:,tr) = psdq;
    end
end


%% -------------------- MOTION ANALYSIS -----------------------------------
% -------------------------------------------------------------------------
%% LOWPASS MOTION
MOTsmooth = MOT;
for tr = 1:Ntr
    MOTsmooth(:,tr) = sgolayfilt(MOT(:,tr), 1, 1*Fs+1);                     % Lowpass motion with Savitzky-Golay filter
end

%% DETECT MOTION BOUTS 
critV1 = 0.02;                                                              % Criterion for movement detection (fixed value)
critV2 = 0.01;                                                              % Criterion of immobility

Move = nan(Ntr,lt);
Move_time = cell(1,Ntr);

for tr = 1:Ntr
    [Move(tr,:),Move_time{tr}] = detect_movement(MOTsmooth(:,tr),Fs,[critV1, critV2]); % DETECT MOTION
end
Move = logical(Move);

%% SPLIT MOTION/IMMOBILITY RATES
Rm = cell(1,Ntr);
Ri = cell(1,Ntr);

for tr = 1:Ntr
    Rm{tr} = [];
    mot = [];
    for m = 1:size(Move_time{tr},1)                                                 % For each motion segment
        [~,k1] = min(abs(bins - Move_time{tr}(m,1)));                               % Find the bin closest to motion initiation
        [~,k2] = min(abs(bins - Move_time{tr}(m,2)));                               % and termination
        k2 = min(k2,lb);                                % In case it goes just above lb
        Rm{tr} = [Rm{tr}, R(:,k1:k2,tr)];                                                  % Concatenate the motion rates
        
        mot = [mot, k1:k2];                                                     % Keep the motion indexes
    end
    
    Ri{tr} = R(:,:,tr);
    Ri{tr}(:,mot) = [];                                                             % Keep concatenated immobility rates
end

%% FREQ-LOCOMOTION CORRELATIONS
% Corr = ones(Nr,lf);
% Pcorr = nan(Nr,lf);
% for r = 1:Nr 
%     for f = 1:lf
%         [Corr(r,f),Pcorr(r,f)] = corr(Ampl(r,:,f)',MOT','type','Pearson');                     % Correlation of theta power and motion over concatenated trials
%     end
% end

% % Get shuffled correlations
% Corr_sh = ones(Nr,reps); 
% for r = 1:Nr
%     lags = 2*rand(1,reps) - 1;                                              % 1000 random numbers in [-1 1] range
%     lags = floor(lags*length(time)/2);                                      % Get random (round) lags up to +- duration/2 of total time
%     for z = 1:reps 
%         Corr_sh(r,z) = corr(Ampl(r,:)',circshift(MOT',lags(z)),'type','Pearson');
%     end
% end

%% MOTION/IMMOBILITY SPECTRA (CONCATENATED TRIALS)
Sm = nan(Nr,8001);
Si = nan(Nr,8001);                                                          % 4001 = length of Freq vector

Moveconcat = Move';
Moveconcat = Moveconcat(:);

despDFFconcat = reshape(permute(despDFF,[2,3,1]),lt*Ntr,[]);            % Concatenate all trials in one column
DFFm = despDFFconcat(Moveconcat',:);                                                         % Keep traces during motion
DFFi = despDFFconcat(~Moveconcat',:);                                                        % and immobility

Fvec = nan;
if sum(Moveconcat) >= 2*Fs         % Change this to 3*Fs for 6_16 PV5                                           % IF THERE ARE AT LEAST 2 sec OF MOTION
    for r = 1:Nr
        [Sm(r,:),Fvec] = pspectrum(DFFm(:,r),Fs,'FrequencyLimits',[0 100],'FrequencyResolution',0.5);      % Compute power spectrum during motion
    end
end

if sum(~Moveconcat) >= 2*Fs                                                   % IF THERE ARE AT LEAST 2 sec OF IMMOBILITY
    for r = 1:Nr
        [Si(r,:),Fvec] = pspectrum(DFFi(:,r),Fs,'FrequencyLimits',[0 100],'FrequencyResolution',0.5);
    end
end

Sm(~isnan(Sm)) = pow2db(Sm(~isnan(Sm)));                                    % Switch to decibel (skip Nan entries where it crashes)
Si(~isnan(Si)) = pow2db(Si(~isnan(Si)));

%% MOTION/IMMOBILITY BAND POWER
BPm = nan(Nr,Ntr,lf);
BPi = nan(Nr,Ntr,lf);

for tr = 1:Ntr
    if sum(Move(tr,:)) >= 1*Fs                                                    % IF THERE ARE AT LEAST 1 sec OF MOTION
        DFFm = despDFF(:,Move(tr,:),tr);                                                         % Keep traces during motion
        for r = 1:Nr
            for f = 1:lf
                BPm(r,tr,f) = bandpower(DFFm(r,:),Fs,freq(f,:));              % Get average theta power during motion
            end
        end
    end
    if sum(~Move(tr,:)) >= 1*Fs                                                   % IF THERE ARE AT LEAST 1 sec OF IMMOBILITY
        DFFi = despDFF(:,~Move(tr,:),tr);                                                        % and immobility
        for r = 1:Nr
            for f = 1:lf
                BPi(r,tr,f) = bandpower(DFFi(r,:),Fs,freq(f,:));
            end
        end
    end
end


%% ORGANIZE OUTPUT
V.DFF = DFF;
V.DespDFF = despDFF;
V.R = R;
V.SP = SP;
V.Template = Template;
V.ROI = ROI;
V.CC = CC;
V.time = time;
V.bins = bins;

B.pars = pars;
B.trials = trials;
B.outcomes = outcomes;
B.progress = progress;
B.licks = licks;

F.PSD = PSD;
F.Phase = Phase;
F.Ampl = Ampl;
F.SP_ph = SP_ph;
F.pvalR = pvalR;
F.pvalV = pvalV;
F.mDir = mDir;
F.lVec = lVec;
F.Skew = Skew;
F.freq = freq;
F.tq = tq;
F.fq = fq;

S.ISI = ISI;
S.PeriSp = PeriSp;
S.Trough = Trough;
S.PTdist = PTdist;
S.Halfwidth = Halfwidth;
S.Points = Points;
S.mR = mR;
S.BI = BI;
S.CSI = CSI;

M.MOT = MOT;
M.Move = Move;
M.Move_time = Move_time;
% M.Corr = Corr;
% M.Pcorr = Pcorr;
M.Sm = Sm;
M.Si = Si;
M.Fvec = Fvec;
M.BPm = BPm;
M.BPi = BPi;
M.Ri = Ri;
M.Rm = Rm;

%% SAVE
datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
save(datafile,'V','B','F','S','M');

