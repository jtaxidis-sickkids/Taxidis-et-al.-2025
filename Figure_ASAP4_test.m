
%% ASAP4 EXAMPLES
% Plot_VOLPY_processing('ASAP4M-B','02_28_2020',1);         % no dip
% % Plot_VOLPY_processing('ASAP4M-B','03_10_2020',5);         % no dip
% % Plot_VOLPY_processing('ASAP4M-B','03_10_2020',7);         % no dip
% % Plot_VOLPY_processing('ASAP4M-B','03_11_2020',1);         % POSITIVE DEFLECTION (but doesnt look like artifact)
% % Plot_VOLPY_processing('ASAP4M-B','03_12_2020',1);         % no dip
% 
% Plot_VOLPY_processing('ASAP4G-A','02_12_2020',1);         % small dip
% Plot_VOLPY_processing('ASAP4G-A','02_12_2020',2);         % shows dip
% % Plot_VOLPY_processing('ASAP4G-A','02_12_2020',4);         % no dip
% 
% Plot_VOLPY_processing('MM6-1','08_20_2020',3);         % shows dip
% Plot_VOLPY_processing('MM6-1','08_20_2020',5);         % shows dip
% % Plot_VOLPY_processing('MM6-2','08_20_2020',1);         % no dip
% % Plot_VOLPY_processing('MM6-3','08_21_2020',2);         % no dip
% % Plot_VOLPY_processing('MM6-3','08_21_2020',5);         % no dip

%% DAYS WITH VALVE IN HAND
%       ID           day        video  maxtrial 
days = {'ASAP4M-B', '02_28_2020', 1,    8;         
        'ASAP4G-A', '02_12_2020', 1,    10;         
        'ASAP4G-A', '02_12_2020', 2,    10;         
           'MM6-1', '08_20_2020', 3,    8;         
           'MM6-1', '08_20_2020', 5,    8};

ld = size(days,1);

%% PLOT EXAMPLE 
Plot_VOLPY_processing('MM6-1','08_20_2020',3);         % shows dip

%% LOAD STANDARD STUFF
asapfile = get_ASAPfile('SlideTest','02_18_2022',4);
[path,videoname] = fileparts(asapfile);
load(fullfile(path,[videoname,'_data.mat']));
time = V.time;
pars = B.pars;
timepoints = [pars.stim1_on, pars.delay_on, pars.stim2_on, pars.stim2_off, ...
    pars.lick_on, pars.vacuum_on,  pars.trial_dur];
lt = length(time);

%% GET MEAN DFF
Fs = 1000;
time = 0 : 1/Fs : (11*Fs-1)/Fs;
MD = nan(ld,lt);
SD = nan(ld,lt);

for d = 1:ld
    [asapfile,~] = get_ASAPfile(days{d,1},days{d,2},days{d,3});
    [path,videoname] = fileparts(asapfile);
    
    datafile = fullfile(path,[videoname,'_data.mat']);
    load(datafile);

    dff = squeeze(V.DFF(1,:,:))';                                           % Keep DFF NOT DESPIKED [Trials x time]
    dff = dff(1:days{d,4},:);                                               % Keep up to max trial
    dff = sgolayfilt(dff, 1, 101, [], 2);                                   % Filter
    MD(d,:) = mean(dff,1);
    SD(d,:) = SEM(dff);
end

%% ASAP4 MEAN DFFs
figure;  hold on;
cols = 0.3*ones(ld,3);
plot_mrate_traces(MD,SD,time,timepoints,cols);


return


%% COMPARE DFF DIP BEFORE AND AFTER BLEACHING
Plot_VOLPY_processing('MM6-1','08_20_2020',3);         %
Plot_VOLPY_processing('MM6-1','08_20_2020',4);         %
figure;
for d = 1:2
    subplot(2,2,d)
    dff = load_dff('MM6-1','08_20_2020',2+d,1);
    mdff = mean(dff,1);
    sdff = SEM(dff);
    plot_mrate_traces(mdff,sdff,time,timepoints,'k')
end
subplot(221)
title('MM6-1, 08_20_2020, before bleaching');
subplot(222)
title('MM6-1, 08_20_2020, after bleaching');



