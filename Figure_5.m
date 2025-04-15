%% SHOW EXAMPLE FIRSTSPIKE 
Plot_VOLPY_processing('PV1','06_16_2021',8)
Plot_VOLPY_processing('PV1','06_14_2021',7)

%% PLOT SPIKING
t_firstspike = 1 + [0.03 0.05];                                             % LOOK FOR FIRST SPIKE AT 30-50ms AFTER ODOR ONSET
t_rebound = 1 + [0.2 0.5];                                                  % LOOK FOR REBOUND AT 200-300ms AFTER ODOR ONSET
[MRpv,~] = Figure_fine_timescale('PV',t_firstspike,t_rebound);
[MRsom,bins] = Figure_fine_timescale('SOM',t_firstspike,t_rebound);

%% PLOT ALL INs SPIKING (REVISIONS)
MR = [cell2mat(MRpv(:,3)); cell2mat(MRsom(:,3))];
figure; hold on
histogram('BinEdges',bins,'BinCounts',mean(MR,1),'EdgeColor','none','FaceColor','k');
fill_plot(bins(1:end-1),mean(MR,1),SEM(MR),'k'); % Plot mean+SEM curves

axis tight;
xlim([0.8 2.2]);
ylabel('All cells all trials');

%% COMPARE FIRST SPIKING AND REBOUND IN PV VS SOM
MR1 = cell2mat(MRpv(:,3));                                                  % Keep all cells over all trials
MR2 = cell2mat(MRsom(:,3));

baseline = find(bins >= 0.5 & bins <= t_firstspike(1));
mr1_b = mean(MR1(:,baseline),2);
mr2_b = mean(MR2(:,baseline),2);

figure;
for i = 1:2                                                                 % For either firstspikes or rebounds
    if i == 1
        onbins = find(bins >= t_firstspike(1) & bins <= t_firstspike(2));
    else
        onbins = find(bins >= t_rebound(1) & bins <= t_rebound(2));
    end
    
    mr1 = MR1(:,onbins);                                                    % KEep mean rates over selsected bins
    mr2 = MR2(:,onbins);
    
    mr1 = mean(mr1,2);                                                      % Get mean per cell over all the selected bins
    mr2 = mean(mr2,2);
    
    mr1 = 100*(mr1 - mr1_b) ./ mr1_b;                                           % Relative increase compared to baseline
    mr2 = 100*(mr2 - mr2_b) ./ mr2_b;
    
    p = ranksum(mr1,mr2);
    
    subplot(1,2,i); hold on
    hw = table(mr1',mr2','VariableNames',{'PV','SOM'});
    violinplot(hw);
    plot_significance(p,1,2,mr1,mr2)
end

subplot(121); 
ylabel('Mean rate at first spike');
subplot(122); 
ylabel('Mean rate at rebound');

%% PLOT FINE-TIMESCALE ODOR RESPONSES IN ON-OFF TRIALS
Figure_fine_timescale_ONOFF('PV');
Figure_fine_timescale_ONOFF('SOM');

%% PLOT WITH VS WITHOUT SPIKING (FIRST SPIKE DOES NOT DRIVE HYPERPOL)
t_firstspike = 1 + [0.03 0.05];                                             % LOOK FOR FIRST SPIKE AT 30-50ms AFTER ODOR ONSET
t_rebound = 1 + [0.2 0.5];                                                  % LOOK FOR REBOUND AT 200-300ms AFTER ODOR ONSET
                                                   
Figure_FirstSpike('PV',t_firstspike)
Figure_FirstSpike('SOM',t_firstspike)   

%% PROGRESS OF FIRST SPIKE AND REBOUND OVER DAYS
Figure_Progress_Spikes('PV',t_firstspike)    
Figure_Progress_Spikes('SOM',t_firstspike)    

Figure_Progress_Spikes('PV',t_rebound)    
Figure_Progress_Spikes('SOM',t_rebound)    
