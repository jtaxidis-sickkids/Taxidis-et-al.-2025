
%% --------- FIGURE 1 -----------------------------------
% -----RECORDINGS AND SPIKING FEATURES ------------------

%% EXAMPLE CELLS 
Plot_VOLPY_processing('PV1','07_01_2021',2)
Plot_VOLPY_processing('PV1','07_06_2021',5)
Plot_VOLPY_processing('PV3','06_15_2021',1)
Plot_VOLPY_processing('PV3','07_05_2021',1)
Plot_VOLPY_processing('PV5','07_06_2021',1)
Plot_VOLPY_processing('PV5','07_06_2021',7)

Plot_VOLPY_processing('MM7-1','09_21_2020',2)
Plot_VOLPY_processing('MM7-2','08_28_2020',1)
Plot_VOLPY_processing('MM7-2','08_28_2020',3)
Plot_VOLPY_processing('ASAP3-1','05_11_2021',5)
Plot_VOLPY_processing('ASAP3-3','04_24_2021',5)
Plot_VOLPY_processing('ASAP3-4','04_13_2022',2)

%% EXAMPLE LONG RECORDINGS
Plot_VOLPY_processing('PV1','06_16_2021',1)
Plot_VOLPY_processing('PV1','07_01_2021',3)
Plot_VOLPY_processing('PV1','07_03_2021',3)

Plot_VOLPY_processing('ASAP3-3','05_11_2021',2)

%% EXAMPLE MULTIPLE CELLS (NOT USED)
Plot_VOLPY_processing('PV5','06_25_2021',4)
Plot_VOLPY_processing('ASAP3-3','05_08_2021',7)

%% COMPARE PV - SOM
Figure_PVvsSOM_spikes;
Figure_PVvsSOM_freq;

%% SPIKE ANALYSIS
[Rm{1},Ri{1}] = Plot_UnitAnalysis_pooled('PV');
[Rm{2},Ri{2}] = Plot_UnitAnalysis_pooled('SOM');

%% COMPARE THETA AND RATES DUIRING LOCOMOTION AND ODOR (WIHTOUT EDGE EFFECTS) (REVISIONS)
Theta_Rates_Locomotion('PV');
Theta_Rates_Locomotion('SOM');

%% COMPUTE RATES DURING REWARD ZONE IMMOBILITY (REVISIONS)
Ri_rew{1} = Reward_Immobility_Rates('PV');
Ri_rew{2} = Reward_Immobility_Rates('SOM');

%% COMPARE PV/SST MOTION & IMMOBILITY RATES (REVISIONS)
[~,pRm] = ttest2(Rm{1},Rm{2});           % 2 sample ttest
[~,pRi] = ttest2(Ri{1},Ri{2});           % 2 sample ttest
[~,pRi_rew] = ttest2(Ri_rew{1},Ri_rew{2});           % 2 sample ttest

figure;
subplot(131);
hw = table(Rm{1}',Rm{2}','VariableNames',{'PV','SST'});
violinplot(hw);
plot_significance(pRm,1,2,Rm{1},Rm{2});
title('Locomotion rates');
subplot(132);
hw = table(Ri{1}',Ri{2}','VariableNames',{'PV','SST'});
violinplot(hw);
plot_significance(pRi,1,2,Ri{1},Ri{2});
title('Immobility rates');
subplot(133);
hw = table(Ri_rew{1}',Ri_rew{2}','VariableNames',{'PV','SST'});
violinplot(hw);
plot_significance(pRi_rew,1,2,Ri_rew{1},Ri_rew{2});
title('Reward Immobility rates');

%% RATES-LOCOMOTION CORRELATION (USE ALSO FOR FIGURE 2)
Speed_Score('PV');
Speed_Score('SOM');

%% PLOT JUST THETA (REVISIONS)
% Plot_theta('PV3','06_14_2021',3); % trial 5
Plot_theta('PV3','06_15_2021',1) % trial 1
Plot_theta('PV5','07_06_2021',7); % trial 3
Plot_theta('MM7-2','08_28_2020',3); % trial 1
Plot_theta('ASAP3-3','05_11_2021',2); % trial 4

%% COUNT MULTIPLE ROI
for INclass = {'PV', 'SOM'}
    load(['PooledData_',cell2mat(INclass),'.mat'],'R');
    
    R = vertcat(R{:});
    roi = cellfun(@(x) size(x,1),R);
    roi = roi > 1;
    sum(roi)
end

%% COUNT AVERAGE TRIALS
for INclass = {'PV', 'SOM'}
    load(['PooledData_',cell2mat(INclass),'.mat'],'R');
    
%     R = R(:,1:2);  % Count only naive sessions
%     R(:,1:2) = []; % Count only trained sessions
    
    R = vertcat(R{:});
    ls = length(R); % Number of sessions
    Ntr = cellfun(@(x) size(x,3),R);
    [mean(Ntr) std(Ntr) ls]
    max(Ntr)
end

%% COUNT CELLS PER MOUSE (REVISIONS)
Days_PV;
ASAP{1} = asap;

Days_SOM;
ASAP{2} = asap;

for ct = 1:2
    asap = ASAP{ct};

    NR = nan(1,length(asap));
 
    for a = 1:length(asap)                                                  % For each animal
        nr = 0;
            
        ID = asap(a).name;
        days = asap(a).sessions(:,1);
        ses = asap(a).sessions(:,3);
        for d = 1:length(days)                                              % For each day
            for s = 1:length(ses{d})                                        % For each video
                asapfile = get_ASAPfile(ID,days{d},ses{d}(s));
                disp(asapfile);
                [path,videoname] = fileparts(asapfile);
                modfile = fullfile(path,[videoname,'_Mcells.mat']);
                load(modfile,'Mcells');
                
                nr = nr + size(Mcells,1);                                  % Count Mcells
            end
        end
        NR(a) = nr;
    end
    [mean(NR) std(NR)]
end

%% COMPARE LED POWER (REVISIONS)
naive_ses_SST = [3, 6, 17, 22, 12];
naive_ses_PV = [12, 8, 9, 15, 7];

LED = readtable('LED power.xlsx');

LED_SOM = [LED.Var1; LED.Var2; LED.Var3; LED.Var4; LED.Var5];
LED_PV  = [LED.Var6; LED.Var7; LED.Var8; LED.Var9; LED.Var10];

[nanmean(LED_PV) nanstd(LED_PV)]
[nanmean(LED_SOM) nanstd(LED_SOM)]
[~,p] = ttest2(LED_PV,LED_SOM)
