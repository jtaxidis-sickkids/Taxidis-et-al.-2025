
%% PLOT EXMAPLE OF FIRING RATE RAMPS
Plot_VOLPY_processing('PV1','06_16_2021',5);
Plot_VOLPY_processing('PV2','06_16_2021',5);
Plot_VOLPY_processing('PV5','06_25_2021',5);
Plot_VOLPY_processing('ASAP3-3','04_24_2021',11);
Plot_VOLPY_processing('ASAP3-3','05_08_2021',5);



%% (PLOT ALL CELLS FROM AN ANIMAL)
Days_PV
a = 5;
for d = 1:size(asap(a).sessions,1)
    for s = asap(a).sessions{d,3}
        Plot_over_trials_rate(asap(a).name,asap(a).sessions{d,1},s,1,-1);
    end
end

%% PLOT NEGATIVE SEQUENCE
[lmc,Nr] = Pool_Neg_sequences('PV')         
[lmc,Nr] = Pool_Neg_sequences('SOM')         

