
%% --------- FIGURE 2 ------------------------------------------
% ----- ODOR-CELLS AND THEIR SPIKING FEATURES ------------------

%% EXAMPLE MCELL TRACES
Plot_VOLPY_processing('PV1','06_16_2021',7);     % USE TRIALS 3,5,7,8
Plot_VOLPY_processing('ASAP3-1','04_21_2021',5); % USE TRIALS 1-4

%% EXAMPLE MCELLS RATES
Plot_over_trials_rate('PV1','06_16_2021',7,1);
Plot_over_trials_rate('PV1','06_16_2021',8,1);
Plot_over_trials_rate('PV5','06_16_2021',3,1);

Plot_over_trials_rate('ASAP3-1','04_21_2021',5,1);
Plot_over_trials_rate('ASAP3-1','05_13_2021',6,1);
% Plot_over_trials_rate('ASAP3-3','05_08_2021',5,1);
Plot_over_trials_rate('ASAP3-4','04_01_2022',2,1);

%% EXAMPLE ODOR-SPECIFIC MCELLS
Plot_over_trials_rate('PV1','07_01_2021',1,1);
Plot_over_trials_rate('PV5','07_06_2021',7,1);

Plot_over_trials_rate('ASAP3-1','05_18_2021',2,1);
Plot_over_trials_rate('ASAP3-3','04_24_2021',15,1);

%% (PLOT ALL CELLS FROM AN ANIMAL)
Days_PV
a = 5;
for d = 1:size(asap(a).sessions,1)
    for s = asap(a).sessions{d,3}
        Plot_over_trials_rate(asap(a).name,asap(a).sessions{d,1},s,1,1);
    end
end

%% PLOT ODOR-SPECIFIC AND NON-ODOR-SPECIFIC SEQUENCE
[lmc,Nr] = Pool_sequences('PV')                                             % (With REVISIONS)
[lmc,Nr] = Pool_sequences('SOM')         

%% COMPARE %MCELLS AND SI
Figure_PVvsSOM_Mcells;                                                      % (With REVISIONS)

%% COMPARE CELLS WITH VS WITHOUT FIELDS (REVISIONS)
Figure_field_vs_nofield('PV');
Figure_field_vs_nofield('SOM');

%% COMPARE MCELL RATES WITH VS WITHOUT LOCOMOTION (REVISIONS)
Figure_field_locomotion('PV');
Figure_field_locomotion('SOM');

%% COMPARE RATES BETWEEN ODORS AND DELAY VS BASELINE
Figure_Rates_stats('PV')
Figure_Rates_stats('SOM')

%% SVM
SVM_Decoding('PV');
SVM_Decoding('SOM');

SVM_Decoding('PV','odorbins');
SVM_Decoding('SOM','odorbins');

%% PLOT BAYESIAN DECODING
Bayesian_Decoding('PV');   
Bayesian_Decoding('SOM');   

%% COMPARE RATES IN ON-OFF TRIALS
Figure_ONOFF('PV')
Figure_ONOFF('SOM')

