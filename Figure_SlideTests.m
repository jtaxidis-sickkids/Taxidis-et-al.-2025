
%% LOAD STANDARD STUFF
asapfile = get_ASAPfile('SlideTest','02_18_2022',4);
[path,videoname] = fileparts(asapfile);
load(fullfile(path,[videoname,'_data.mat']));
time = V.time;
pars = B.pars;
timepoints = [pars.stim1_on, pars.delay_on, pars.stim2_on, pars.stim2_off, ...
    pars.lick_on, pars.vacuum_on,  pars.trial_dur];

%% PLOT BEADS
Plot_VOLPY_processing('SlideTest','02_22_2022',1);         % Beads1 = Airflow on the slide
Plot_VOLPY_processing('SlideTest','02_22_2022',2);         % Beads2 = Airflow away from the slide

%% PLOT SLIDE - DAY 1
Plot_VOLPY_processing('SlideTest','02_18_2022',2);         % Slide - 6 trials without foam, 6 trials with
Plot_VOLPY_processing('SlideTest','02_18_2022',3);         % Slide - 6 trials with foam and airflow on slide, 
                                                           %         6 trials with airflow away
                                                           %         6 trials with airflow off      
Plot_VOLPY_processing('SlideTest','02_18_2022',4);         % Slide - 8 trials with foam and airflow on slide, 
                                                           %         5 trials with airflow away
                                                           %         5 trials with olfactometer off      

%% PLOT SLIDE - DAY 2
Plot_VOLPY_processing('SlideTest','02_22_2022',3);         % Slide1 = Airflow on the slide
Plot_VOLPY_processing('SlideTest','02_22_2022',4);         % Slide2 = Airflow away from the slide
 
%% PLOT SLIDE - DAY 3
Plot_VOLPY_processing('SlideTest','02_24_2022',1);         % Slide - 12 trials with slide stuck to headbar holder with tape, 
                                                           % foam, airflow away from slide  

%% PLOT BEADS MEAN DFF
bead = 1;

figure; 
subplot(421)
dff = load_dff('SlideTest','02_22_2022',1,bead);
mdff = mean(dff,1);
sdff = SEM(dff);
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
title('Beads-airflow on slide')

subplot(423)
dff = load_dff('SlideTest','02_22_2022',2,bead);
mdff = mean(dff,1);
sdff = SEM(dff);
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
title('Beads-airflow away from slide')

%% PLOT SLIDES - DAY 2 MAIN DFF
roi = 1;

subplot(425)
dff = load_dff('SlideTest','02_22_2022',3,roi);
mdff = mean(dff,1);
sdff = SEM(dff);
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
title('Slide-airflow away from slide')

subplot(427)
dff = load_dff('SlideTest','02_22_2022',4,roi);
mdff = mean(dff(1:3,:),1);
sdff = SEM(dff(1:3,:));
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
mdff = mean(dff(4:6,:),1);
sdff = SEM(dff(4:6,:));
plot_mrate_traces(mdff+0.5,sdff,time,timepoints,'b')
ylim([-0.1 0.6]);
title('Slide-airflow away, with/without airvalve sitting on rig')

%% PLOT SLIDE - DAY 1 AND DAY 3 MEAN DFF
roi = 1;
figure
subplot(422)
dff = load_dff('SlideTest','02_18_2022',2,roi);
mdff = mean(dff(1:6,:),1);
sdff = SEM(dff(1:6,:));
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
mdff = mean(dff(7:12,:),1);
sdff = SEM(dff(7:12,:));
plot_mrate_traces(mdff+0.5,sdff,time,timepoints,'b')
ylim([-0.4 0.8]);
title('Slide-without foam / with foam')

subplot(424)
dff = load_dff('SlideTest','02_18_2022',3,roi);
mdff = mean(dff(1:6,:),1);
sdff = SEM(dff(1:6,:));
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
mdff = mean(dff(7:12,:),1);
sdff = SEM(dff(7:12,:));
plot_mrate_traces(mdff+0.5,sdff,time,timepoints,'b')
mdff = mean(dff(13:18,:),1);
sdff = SEM(dff(13:18,:));
plot_mrate_traces(mdff+1,sdff,time,timepoints,'r')
ylim([-0.4 1.2]);
title('Airflow on slide / aiflow away / no airflow')

subplot(426)
dff = load_dff('SlideTest','02_18_2022',4,roi);
mdff = mean(dff(1:8,:),1);
sdff = SEM(dff(1:8,:));
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
mdff = mean(dff(9:13,:),1);
sdff = SEM(dff(9:13,:));
plot_mrate_traces(mdff+0.5,sdff,time,timepoints,'b')
mdff = mean(dff(14:18,:),1);
sdff = SEM(dff(14:18,:));
plot_mrate_traces(mdff+1,sdff,time,timepoints,'r')
ylim([-0.4 1.2]);
title('Airflow on slide / aiflow away / no airflow')

subplot(428)
dff = load_dff('SlideTest','02_24_2022',1,roi);
mdff = mean(dff,1);
sdff = SEM(dff);
plot_mrate_traces(mdff,sdff,time,timepoints,'k')
title('Slide stuck with tape, airflow away')
