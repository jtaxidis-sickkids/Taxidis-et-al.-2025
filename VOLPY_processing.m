function VOLPY_processing(ID,day,session)

%% LOAD VOLPY DATA
[asapfile,asapfiles,video] = get_ASAPfile(ID,day,session);
disp(asapfile);
[edrfile,matfile] = get_Behaviorfiles(ID,day,session,asapfile,asapfiles);   % FIND THE CORRECT BEHAVIOR FILES

[path,videoname] = fileparts(asapfile);
VOLPYfile = fullfile(path,'Volpy_data',[videoname,'_volpy_data.mat']);
load(VOLPYfile,'DFF','ROIs','spikes','weights','mean_im');

[~,Nfr] = size(DFF);

if Nfr == 88001
    DFF(:,end) = [];
end

if Nfr == 220724    % For PV4 - 06_25_21 - t4
    DFF(:,end-723:end) = [];
end

if strcmp(ID,'ASAP3-3') && strcmp(day,'05_08_2021') && session == 1
    ltr = 11e3; 
    k = [0 0 10700 10920 10810 10870 10870 10880 10840 10700 10710];
    for tr = 3:11
        DFF = [DFF(1 : ((tr-1)*ltr + k(tr))), zeros(1,ltr-k(tr)), DFF(((tr-1)*ltr + k(tr)+1) : end)];
        
        z = find(spikes > ((tr-1)*ltr + k(tr)), 1, 'first');
        spikes(z:end) = spikes(z:end) + ltr-k(tr);
    end
end
[Nr,Nfr] = size(DFF);

%% SET TRIAL DURATION
trial_dur = 11; 
if contains(videoname,'long')
    trial_dur = 16;
end
Fs_asap = 1e3;
ls = trial_dur*Fs_asap;                         % Length of section (trial)
ASAPtime = (0:ls-1)/Fs_asap;

Ntr = ceil(Nfr/ls);                             % Number of trials
disp(['Number of trials = ',num2str(Ntr)]);

%% RENAME AND RESHAPE DATA
if iscell(mean_im)
    mean_im = mean_im{1};
    mean_im = reshape(mean_im,[1,size(mean_im)]);
end
Template = squeeze(mean_im(1,:,:));

CC = cell(Nr,1);
for r = 1:Nr
    structBoundaries = bwboundaries(squeeze(ROIs(r,:,:)));                  % Get coordinates of the boundary of the freehand drawn region.
    k = structBoundaries{1};                                                % Get n by 2 array of x,y coordinates.
    CC{r} = flipud(k');
end    

SP = cell(Nr,Ntr);
if ~iscell(spikes)
    sptemp = spikes;
    spikes = {};
    for r = 1:Nr
        spikes{r} = sptemp(r,:);
    end
end

totspikes = 0;
for r = 1:Nr
    for tr = 1:Ntr
        k = spikes{r} > (tr-1)*ls & spikes{r} <= tr*ls;        
        sp_tr = spikes{r}(k) - (tr-1)*ls;
        SP{r,tr} = ASAPtime(sp_tr);
        totspikes = totspikes + length(SP{r,tr});
    end
    SP{r,1}(SP{r,1} < 0.2) = [];                                            % REMOVE SPIKES BEFORE THE FIRST 200ms FROM FIRST TRIAL
end
totspikes

if Nfr < Ntr*ls
    DFF = cat(2,DFF,nan(Nr, Ntr*ls - Nfr));
end
DFF = reshape(DFF,[Nr,ls,Ntr]);
DFF = DFF * 100;                                                            % Turn to percentage

%% COMPUTE FIRING RATES
binl = 0.1;                                                                 % Firing rate time bin 100ms
bins = 0 : binl : trial_dur;

R = nan(Nr,length(bins)-1,Ntr);
for r = 1:Nr
    for tr = 1:Ntr
        R(r,:,tr) = histcounts(SP{r,tr},bins) / binl;                       % Compute firing rate
%         R(r,:,tr) = smooth(R(r,:,tr),7);                                  % And smooth it
    end
end

bins = bins + binl/2;
bins(end) = [];

%% GET BEHAVIOR -----------------------------------
F_edr = 1e3;
EDR_trial_dur = trial_dur*F_edr;    % EDR-trial duration

% LOAD EDR
edr = import_edr(edrfile);                                                  % Import the EDR data
mot = edr(:,3);
od1 = edr(:,4);
od2 = edr(:,5);

mot = (mot - mode(mot))/1000;                                               % Turn to V
od1 = (od1 - mode(od1))/1000/5;                                             % Turn to 0-1 (1 = 5000mV)
od2 = (od2 - mode(od2))/1000/5;
mot = abs(mot);                                                             % Make motion be positive

% GET MATLAB PARAMETERS AND DATA
load(matfile);                                                              % Load behavior file
pars = p;                                                                   % Keep all parameters and inputs
Ntr = inputs.Ntr;
combos = progress(:,1);
trials = floor(combos/10);
outcomes = progress(:,2);

% ARRANGE MOTION INTO TRIALS
bothodors = sum(od1+od2,2);                                                 % Add EDR odor signals to get all odor stimuli
cross = (bothodors > 0.1);                                                  % Find when the odors are ON
tr_on = find((cross - circshift(cross,1)) == 1);                            % Find when they turn ON

MOT = nan(EDR_trial_dur,Ntr);
for tr = 1:Ntr
    MOT(:,tr) = mot(tr_on(tr) + (0:EDR_trial_dur-1));
end

%% ORGANIZE FINAL STRUCTURES
V.Template = Template;
V.weights = weights;   
V.ROI = ROIs;
V.CC = CC;
V.DFF = DFF;
V.time = ASAPtime;
V.SP = SP;
V.R = R;
V.bins = bins;

B.MOT = MOT;
B.trials = trials;
B.outcomes = outcomes;
B.combos = combos;
B.pars = pars;
B.progress = progress;
B.Data = Data;

%% STORE DATA 
outfile = fullfile(path,[videoname,'_data.mat']);
save(outfile,'V','B'); 

%% PLOT PROCESSING
Plot_VOLPY_processing(ID,day,session);

