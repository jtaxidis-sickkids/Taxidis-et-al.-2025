function [DFF,SP,R,MOT,trials,outcomes,progress,licks,roi] = load_and_pool(ID,day,session)

%% FIND ID AND DAY INDEXES 
if contains(ID,'PV')
    Days_PV;
else
    Days_SOM;
end

allids = {asap.name};
id = find(ismember(allids,ID));
day_indx = find(ismember(asap(id).sessions(:,1),day));
ses_indx = find(asap(id).sessions{day_indx,3} == session);

pvideos = asap(id).sessions{day_indx,2};

maxtrials = asap(id).sessions{day_indx,4};
maxtrial = maxtrials(ses_indx);

rois = asap(id).sessions{day_indx,5};
if isempty(rois)
    rois = cell(size(maxtrials));
end
roi = rois{ses_indx};

%% KEEP CORRECT VIDEOS
if isempty(pvideos)
    pvideos = session;
else
    k = find(pvideos(:,1) == session);    
    if isempty(k)
        pvideos = session;
    else
        pvideos = pvideos(k,:);
        pvideos = pvideos(1) : pvideos(2);
    end
end

if strcmp(ID,'PV4') && strcmp(day,'06_29_2021') && session == 1
   pvideos = [1 4]; 
end

Nv = length(pvideos);

if Nv > 1
    disp('POOLING VIDEOS:');
end

%% POOL VIDEOS
DFF = [];
SP = {};
R = [];
MOT = [];
trials = [];
outcomes = [];
progress = [];
licks = [];

for v = 1:Nv
    asapfile = get_ASAPfile(ID,day,pvideos(v));
    disp(asapfile);

    [path,videoname] = fileparts(asapfile);
    Datafile = fullfile(path,[videoname,'_data.mat']);
    
    load(Datafile);
    
    DFF = cat(3,DFF,V.DFF);
    SP = cat(2,SP,V.SP);
    R = cat(3,R,V.R); 
        
    MOT = cat(2,MOT,B.MOT);
    trials = cat(1,trials,B.trials);                                                
    outcomes = cat(1,outcomes,B.outcomes);
    progress = cat(1,progress,B.progress);
    licks = cat(2,licks,squeeze(B.Data(:,2,:)));
    % pars = ....
end

if Nv > 1
    disp('-------');
end

%% KEEP ONLY SELECTED ROI
if isempty(roi)
    roi = 1:size(DFF,1);
end
DFF = DFF(roi,:,:);
SP = SP(roi,:);
R = R(roi,:,:);

%% KEEP UP TO MAX TRIAL
DFF = DFF(:,:,1:maxtrial);
SP = SP(:,1:maxtrial);
R = R(:,:,1:maxtrial);
MOT = MOT(:,1:maxtrial);
trials = trials(1:maxtrial);
outcomes = outcomes(1:maxtrial);
progress = progress(1:maxtrial,:);
licks = licks(:,1:maxtrial);

