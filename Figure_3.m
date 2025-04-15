%% --------- FIGURE 3 part 1 -----------------------------------
% ----- STABILITY ACROSS DAYS -------------------------

%% PLOT EXAMPLE STABLE CELL TRACES
Plot_VOLPY_processing('PV1','06_14_2021',7);
Plot_VOLPY_processing('PV1','06_16_2021',7);
Plot_VOLPY_processing('PV1','06_27_2021',4);

Plot_VOLPY_processing('ASAP3-4','04_01_2022',3);
Plot_VOLPY_processing('ASAP3-4','04_13_2022',3);

%% PLOT EXAMPLE STABLE CELL RATES
Figure_StableCells('SOM');               % ASAP3-3 GROUP 2 - RELATIVELY STABLE FIELD FOR 2 DAYS
                                         % ASAP3-3 GROUP 4 - EXAMPLE UNSTABLE FIELD FOR 2 DAYS
                                         % ASAP3-4 GROUP 2 - STABLE FIELD PRE-POST TRAINING
                                         
Figure_StableCells('PV');                % PV1 - GROUP 5 - EXAMPLE OF CELL RECORDED FOR 5 CONSEC. DAYS - ALTERNATING FIELD
                                         % PV1 - GROUP 1 - SUPER STABLE FIELD PRE-POST TRAINING
                                         % PV1 - GROUP 7 - SUPER STABLE FIELD POST TRAINING
                                         % PV1 - GROUP 4 - RELATIVELY STABLE FIELD
                                         % PV1 - GROUP 2 - RELATIVELY STABLE FIELD PRE-POST TRAINING
                                
                                         % PV5 - GROUP 2 - SUPER STABLE FIELD PRE-POST TRAINING
                                         % PV5 - GROUP 3 - EXAMPLE OF EMERGING FIELD           
%% PLOT POOLED STABLE CELLS
Figure_StableCells_Pooled('PV','pre_post');   
Figure_StableCells_Pooled('PV','trained');   
% Figure_StableCells_Pooled('PV','all');   

Figure_StableCells_Pooled('SOM','pre_post');   
Figure_StableCells_Pooled('SOM','trained');   
% Figure_StableCells_Pooled('SOM','all');   

%% PLOT EVOLUTION OF SI AND CELLTYPE ACROSS DAYS
[PVSTABILITY,PVINFLOW,PVOUTFLOW] = Figure_StableCells_vs_Days('PV');
[SOMSTABILITY,SOMINFLOW,SOMOUTFLOW] = Figure_StableCells_vs_Days('SOM');

pstab = ranksum(PVSTABILITY,SOMSTABILITY)

%% ANOVA
y = [PVSTABILITY';PVINFLOW';PVOUTFLOW';SOMSTABILITY';SOMINFLOW';SOMOUTFLOW'];
g1 = [repmat('PVstable',size(PVSTABILITY'));repmat('PVinflow',size(PVINFLOW'));...
    repmat('PVoutflo',size(PVOUTFLOW'));repmat('SMstable',size(SOMSTABILITY'));...
    repmat('SMinflow',size(SOMINFLOW'));repmat('SMoutflo',size(SOMOUTFLOW'))];
[p,~,stats] = anova1(y,g1)

%% PLOT FIRING RATES, THETA AMPL AND VECTOR LENGTH VS DISTANCE BETWEEN DAYS AND IN PRE-POST/X-X+1
[R1pv,R2pv,A1pv,A2pv,L1pv,L2pv] = Plot_Ampl_vs_Days('PV');
[R1som,R2som,A1som,A2som,L1som,L2som] = Plot_Ampl_vs_Days('SOM');

R1 = [R1pv;R1som];
R2 = [R2pv;R2som];
A1 = [A1pv;A1som];
A2 = [A2pv;A2som];   
L1 = [L1pv;L1som];
L2 = [L2pv;L2som];
    
[~,PR_prepost] = ttest(R2(:,1),R2(:,2)) 
[~,PR_post] = ttest(R1(:,1),R1(:,2))                                          
[~,PA_prepost] = ttest(A2(:,1),A2(:,2)) 
[~,PA_post] = ttest(A1(:,1),A1(:,2))                                           
[~,PL_prepost] = ttest(L2(:,1),L2(:,2)) 
[~,PL_post] = ttest(L1(:,1),L1(:,2))                                            

%% PLOT FIRING RATE CORRELATIONS VS DISTANCE BETWEEN DAYS
Plot_RateCorr_vs_Days('PV');
Plot_RateCorr_vs_Days('SOM');



%% --------- FIGURE 3 part 2 -----------------------------------
% ----- PROGRESS ACROSS DAYS --------------

%% PLOT PERFORMANCE, LOCOMOTION AND MCELL FEATURES ACROSS DAYS AND LEARNING
Perf_pv = Figure_Progress('PV','all');
Figure_Progress('PV','mcells');
Figure_Progress('PV','odorcells');
Figure_Progress('PV','nonodorcells');
Figure_Progress('PV','nonmcells');

%% PLOT PERFORMANCE, LOCOMOTION AND MCELL FEATURES ACROSS DAYS AND LEARNING
Perf_sst = Figure_Progress('SOM','all');
Figure_Progress('SOM','mcells');
Figure_Progress('SOM','odorcells');
Figure_Progress('SOM','nonodorcells');
Figure_Progress('SOM','nonmcells');

%% COMPARE PVvsSST PERFORMANCES (REVISIONS)
P1 = cell2mat(Perf_pv(:));
P2 = cell2mat(Perf_sst(:));

[mean(P1) std(P1)/sqrt(length(P1))]
[mean(P2) std(P2)/sqrt(length(P2))]

[~,p] = ttest2(P1,P2)
figure; hold on
hw = table(P1',P2','VariableNames',{'PV','SOM'});
violinplot(hw);
plot_significance(p,1,2,P1,P2)
