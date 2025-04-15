

%% PLOT POOLED FREQ AND MOTION DATA 
Plot_UnitAnalysis_pooled('SOM');

%% PLOT EACH CELL'S ACTIVITY
Days_PV;

for a = 3%:length(asap)
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);  
    for d = 1:length(days)
        for s = 1:length(ses{d})
%             Plot_over_trials_rate(ID,days{d},ses{d}(s),1);
            Plot_over_trials_DFF(ID,days{d},ses{d}(s),1);
%             Plot_over_trials_ampl(ID,days{d},ses{d}(s),1);
%             Plot_over_trials_phase(ID,days{d},ses{d}(s),1);
        end
    end
end

%% PLOT ODOR-SPECIFIC AND NON-ODOR-SPECIFIC SEQUENCE
[lmc,Nr] = Pool_sequences('PV')         
                                       
%% PLOT BANDPOWER OF MCELLS
Pool_sequences_DFF('PV');

%% PLOT BANDPOWER OF MCELLS
Pool_sequences_Bandpower('SOM');

%% PLOT PHASES OF MCELLS
Pool_sequences_Phases('PV');  % PHASE RESETTING. RESETTING DROPS AFTER LEARNIGN (OR WITH DAYS)

%% 
Figure_Progress_DFF;   % SLOW DECREASE OF DFF DROP OVER DAYS ONLY FOR SOM
Figure_Progress;        % Different trends in PV vs SOM