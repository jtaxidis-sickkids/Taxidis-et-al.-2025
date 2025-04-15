%% --------- FIGURE 4 -----------------------------------
% ----- DFF FEATURES & BANDPOWER -------------------------

%% SHOW EXAMPLES OF THE DFF DIP
Plot_VOLPY_processing('PV1','06_14_2021',7)
Plot_VOLPY_processing('PV3','06_14_2021',3)
Plot_VOLPY_processing('PV5','07_06_2021',7)

Plot_VOLPY_processing('MM7-1','09_04_2020',3)
% Plot_VOLPY_processing('ASAP3-3','04_24_2021',5)
Plot_VOLPY_processing('ASAP3-3','05_11_2021',2)
Plot_VOLPY_processing('ASAP3-4','04_13_2022',8)

%% PLOT THE DESPIKED DFF OF SAME CELLS ACROSS TRIALS (SHOW DIP OF AVERAGE DFF)
Plot_over_trials_DFF('PV1','06_14_2021',7,1);
Plot_over_trials_DFF('PV3','06_14_2021',3,-1);
Plot_over_trials_DFF('PV5','07_06_2021',7,1);

Plot_over_trials_DFF('MM7-1','09_04_2020',3,-1);
Plot_over_trials_DFF('ASAP3-3','05_11_2021',2,-1);
Plot_over_trials_DFF('ASAP3-4','04_13_2022',8,1);

%% PLOT SPECTROGRAMS AND FILTERED DFF OF SAME CELLS (SHOW DELTA AND THETA)
Plot_UnitAnalysis('PV1','06_14_2021',7);
Plot_UnitAnalysis('PV3','06_14_2021',3);
Plot_UnitAnalysis('PV5','07_06_2021',7);

Plot_UnitAnalysis('MM7-1','09_04_2020',3);
Plot_UnitAnalysis('ASAP3-3','05_11_2021',2);
Plot_UnitAnalysis('ASAP3-4','04_13_2022',8);

%% DFF MAX AND MIN STATS AND DISTRIBUTIONS
[MinPV,MaxPV,MinOnPV,MinDurPV,MinTPV] = Figure_DFF_stats('PV');
[MinSOM,MaxSOM,MinOnSOM,MinDurSOM,MinTSOM] = Figure_DFF_stats('SOM'); 

%% COMPARE STATS
p = nan(5,1);
p(1) = ranksum(cellfun(@nanmean, MinPV),cellfun(@nanmean, MinSOM));
p(2) = ranksum(cellfun(@nanmean, MaxPV),cellfun(@nanmean, MaxSOM));
p(3) = ranksum(cell2mat(MinOnPV),cell2mat(MinOnSOM));
p(4) = ranksum(cell2mat(MinTPV),cell2mat(MinTSOM));
p(5) = ranksum(cell2mat(MinDurPV),cell2mat(MinDurSOM));

for i = 1:length(MinTPV)                                                    % For each cell     
    sigDipPV{i} = ~isnan(MinTPV{i});                                            % keep "SIGNIFICANT" DIPS indexes
end
for i = 1:length(MinTSOM)                                                    % For each cell
    sigDipSOM{i} = ~isnan(MinTSOM{i});                                            % keep "SIGNIFICANT" DIPS indexes
end
sigdippv = cellfun(@nansum,sigDipPV) ./  cellfun(@length,sigDipPV) * 100;            % Ratio of significant dips per cell
sigdipsom = cellfun(@nansum,sigDipSOM) ./  cellfun(@length,sigDipSOM) * 100;            % Ratio of significant dips per cell
p(6) = ranksum(sigdippv,sigdipsom);

[~,p] = fdr(p);
pmin = p(1)
pmax = p(2)
pon = p(3)
pmint = p(4) 
pdur = p(5) 
psig = p(6)

%% PLOT POOLED AVERAGE DFFs AND FILTERED DFFs OF CELL GROUPS
[MD{1},FMD{1},Fp{1}] = Pool_sequences_DFF('PV');                       
[MD{2},FMD{2},Fp{2}] = Pool_sequences_DFF('SOM');

%% PLOT AVERAGE DFF OF ALL CELLS IN ALL TRIALS (REVISIONS)
load('PooledData_PV.mat','timepoints','time');
MDall = [cell2mat(MD{1}(:,3)); cell2mat(MD{1}(:,3))];
figure;
plot_mrate_traces(nanmean(MDall,1),SEM(MDall),time,timepoints,'k');       % Plot mean across cells    
axis tight;
xlim([0.8 2.2]);
ylabel('All cells all trials');

%% PLOT MEAN DFF OF ODOR AND TIME CELLS SEPARATELY (REVISIONS - FIGURE S2)
load('PooledData_PV.mat','timepoints','time');
ontime = find(time >= timepoints(1) & time < timepoints(2));                % COMPARE DFFS OVER 1ST ODOR 
tontime = time(ontime);

ylabels = {'Odor-spec';'Non-spec';'Non-field'};
tits = {'Preferred';'All trials';'All trials'};
    
figure('Name','Average DFFs');
for ct = 1:2
    for i = 1:3                                                             % For odor-specific / non specific / no-field cells
        subplot(3,2, (i-1)*2 + ct); hold on;
        if i == 1
            A = [MD{ct}{1,1}; MD{ct}{2,2}];                                 % KEEP PREFERRED TRIALS OF ODOR SPECIFIC CELLS
            F = [Fp{ct}{1}; Fp{ct}{2}];
        elseif i == 2
            A = MD{ct}{3,3};                                                % KEEP ALL TRIALS OF NON-SPECIFIC CELLS
            F = Fp{ct}{3};
        else
            A = MD{ct}{4,3};                                                % KEEP ALL TRIALS OF NO-FIELD CELLS
            F = Fp{ct}{4};
        end
    
        od = (F <= 2);                                                      % Find odor-fields

        fill_plot(time,nanmean(A(od,:),1),SEM(A(od,:)),'b');
        fill_plot(time,nanmean(A(~od,:),1),SEM(A(~od,:)),'r');
        xlim([0, timepoints(5)]);
        ylabel(ylabels{i});
        title(tits{i});
        
        % SIGNIFICANCE DURING ODOR
        p = nan(1,length(ontime));
        for t = 1:length(ontime)
            p(t) = ranksum(A(od,t),A(~od,t));
        end
        [~,p] = fdr(p);
        plot(tontime(p<0.05),1.1*max(nanmean(A(od,:),1))*ones(sum(p<0.05),1),'*k')
    end
end
%     for tt1 = 1:4                                                  % For each trial type
%         for tt2 = 1:3                                              % For each trial type (again)
%             subplot(4,3,(tt1-1)*3 + tt2); hold on;
%             od = Fp{ct}{tt1} <= 2;
%             fill_plot(time,nanmean(MD{ct}{tt1,tt2}(od,:),1),SEM(MD{ct}{tt1,tt2}(od,:)),'b');
%             fill_plot(time,nanmean(MD{ct}{tt1,tt2}(~od,:),1),SEM(MD{ct}{tt1,tt2}(~od,:)),'r');
%             xlim([0, timepoints(5)]);
%             ylabel(ylabels{tt1});
%             if tt1 == 1, title(tits{tt2}); end
%             
%             % SIGNIFICANCE DURING ODOR + 1 sec DELAY ONLY
%             p = nan(1,length(ontime));
%             for i = 1:length(ontime)
%                 p(i) = ranksum(MD{ct}{tt1,tt2}(od,i),MD{ct}{tt1,tt2}(~od,i));
%             end
%             [~,p] = fdr(p);
%             plot(tontime(p<0.05),1.1*max(nanmean(MD{ct}{tt1,tt2}(od,:),1))*ones(sum(p<0.05),1),'*k')
%         end
%     end

%% PLOT MEAN DFF AND FILTERED DFF OF FIELD CELLS AROUND THEIR FIELD (REVISIONS - FIGURE S2)
window = -1000 : 1000;
ylabels = {'OdorA';'OdorB';'Non-spec';'Non-field'};
tits = {'OdorA','OdorB','All trials'};

for ct = 1:2
    periMD = cell(4,3);
    periFMD = cell(4,3);
    for tt1 = 1:4
        for tt2 = 1:3
            for r = 1:length(Fp{ct}{tt1})
                periMD{tt1,tt2}(r,:)   =  MD{ct}{tt1,tt2}(r,round(Fp{ct}{tt1}(r) * 1000) + window);
                periFMD{tt1,tt2}(r,:,:) = FMD{ct}{tt1,tt2}(r,round(Fp{ct}{tt1}(r) * 1000) + window,:);
            end
        end
    end
    
    figure('Name','Average peri-field DFFs');
%     for tt1 = 1:4                                                  
%         for tt2 = 1:3                                             
%             subplot(4,3,(tt1-1)*3 + tt2); hold on;
%             for r = 1:length(Fp{ct}{tt1})
%                 plot(window/1000,periMD{tt1,tt2}(r,:),'Color',[0.8,0.8,0.8])
%             end
%             fill_plot(window/1000,nanmean(periMD{tt1,tt2},1),SEM(periMD{tt1,tt2}),'k')
%         end
%     end
                                              
   A = [periMD{1,1};periMD{2,2};periMD{3,3}];       % KEEP PREFERRED TRIALS OF ODOR-SPEC AND ALL TRIALS OF NON-SPEC
   B = [periFMD{1,1};periFMD{2,2};periFMD{3,3}];
   
   F = [Fp{ct}{1}; Fp{ct}{2}; Fp{ct}{3}];
   od = find(F <= 2);                                                      % Find odor-fields
   de =  find(F > 2);
   
   subplot(2,2,1); hold on;
   for r = 1:length(od)
       plot(window/1000,A(od(r),:),'Color',[0.8,0.8,0.8])
   end
   fill_plot(window/1000,nanmean(A(od,:),1),SEM(A(od,:)),'k')
   
   subplot(2,2,3); hold on;
   for r = 1:length(de)
       plot(window/1000,A(de(r),:),'Color',[0.8,0.8,0.8])
   end
   fill_plot(window/1000,nanmean(A(de,:),1),SEM(A(de,:)),'k')
   
        
   f = 2
   subplot(2,2,2); hold on;
   for r = 1:length(od)
       plot(window/1000,B(od(r),:,f),'Color',[0.8,0.8,0.8])
   end
   fill_plot(window/1000,nanmean(B(od,:,f),1),SEM(B(od,:,f)),'k')
   
   subplot(2,2,4); hold on;
   for r = 1:length(de)
       plot(window/1000,B(de(r),:,f),'Color',[0.8,0.8,0.8])
   end
   fill_plot(window/1000,nanmean(B(de,:,f),1),SEM(B(de,:,f)),'k')
   
end

%% PLOT AVERAGE SPECTROGRAM (REVISIONS)
Plot_av_spectrogram('PV');
Plot_av_spectrogram('SOM');

%% PLOT BANDPOWER OF SEQUENCES
Pool_sequences_ampl('PV')         
Pool_sequences_ampl('SOM')  

%% PLOT POOLED PHASES OF SEQUENCES AND PHASE VARIANCE
Pool_sequences_Phases('PV')         
Pool_sequences_Phases('SOM')         

%% PLOT DFF DIP ACROSS DAYS
Figure_Progress_DFF('PV')
Figure_Progress_DFF('SOM')

%% PLOT CORRELATION OF HYPERPOLARIZATIONS WITH PERFORMANCE (REVISIONS)
Figure_DFF_Perf('PV')
Figure_DFF_Perf('SOM')

%% PLOT CORRELATION OF HYPERPOLARIZATIONS WITH MOTION/RATES/THETA (REVISIONS)
Figure_DFF_Mot_R_Theta('PV')
Figure_DFF_Mot_R_Theta('SOM')

%% CONTROL TESTS FOR DFF DIP!!!!!!!!
Figure_SlideTests
Figure_ASAP4_test
Figure_ValveInHand 
Figure_ONOFF('PV')          % USE THE DFF PART
Figure_ONOFF('SOM')         % USE THE DFF PART

% show if dip decreases over trials together with SNR

%% PLOT BANDPOWER OF INDIVIDUAL CELLS (SAME AS FOR RATE)
% Plot_over_trials_ampl('PV1','06_14_2021',7,1);
% Plot_over_trials_ampl('PV3','06_14_2021',3,-1);
% Plot_over_trials_ampl('PV5','07_06_2021',7,1);
% 
% Plot_over_trials_ampl('MM7-1','09_04_2020',3,-1);
% Plot_over_trials_ampl('ASAP3-3','05_11_2021',2,-1);
% Plot_over_trials_ampl('ASAP3-4','04_13_2022',8,1);


%% PLOT PHASES AND PHASE LOCKING OF INDIVIDUAL CELLS (SAME AS FOR RATE)
% Plot_over_trials_phase('PV1','06_16_2021',7,1);
% Plot_over_trials_phase('PV1','06_16_2021',8,1);
% Plot_over_trials_phase('PV3','06_29_2021',5,1);
% 
% 
% Plot_over_trials_phase('PV1','06_14_2021',7,1)
% Plot_over_trials_phase('PV3','06_14_2021',3,-1)
% Plot_over_trials_phase('PV5','07_06_2021',7,1)
% 
% Plot_over_trials_phase('MM7-1','09_04_2020',3,1)
% Plot_over_trials_phase('ASAP3-3','05_11_2021',2,1)
% Plot_over_trials_phase('ASAP3-4','04_13_2022',8,1)



