function [R1,R2,A1,A2,L1,L2] = Plot_Ampl_vs_Days(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV
elseif strcmp(INclass,'SOM')
    Days_SOM
end

Na = length(asap);                                                          % total number of animals

%% POOL RATES OVER TWO CONSECUTIVE IMAGING SESSIONS OR DO PRE-POST TRAINING COMPARISON
MR = cell(Na,2);
AMP = cell(Na,2);
LV = cell(Na,2);
Ddist = cell(Na,2);

for a = 1:Na                                                                % For each animal
    ID = asap(a).name;
    ls = length(asap(a).samecells);
    
    for s = 1:ls                                                            % For each group of registered cells
        reg = asap(a).samecells{s};                                         % Keep the info of all registered cells
        days = reg(:,1);                                                    % Keep the day indexes
        
        Days = cell(1,2);
        Days{1} = find(days > 2);                                           % Turn to row vector
        last_pre = find(days <= 2,1,'last');                                % Find the last naive imaging day of the cell (keep index in list of days)
        first_post = find(days > 2,1,'first');                              % and the first trained day
        if ~isempty(last_pre) && ~isempty(first_post)                       % If the cell was imaged before and after training
            Days{2} = [last_pre; first_post];                               % Keep that last naive day
        end

        for i = 1:2
            for d1 = 1 : length(Days{i})-1                                  % For each selected registered day (except last)
                % FIND CORRECT FILE TO LOAD
                D1 = Days{i}(d1);
                Sess1 = reg(D1,1);                                          % Actual imaging session number
                day = asap(a).sessions{Sess1,1};                            % Actual date
                sess = reg(D1,2);                                           % Actual video index
                roi = reg(D1,3);                                            % ROI index
                
                % LOAD CORRECT FILE (CANNOT REPLACE WITH POOLED DATA DUE TO SESSION NUMBERING IN DAYS_PV/SOM)
                [asapfile,~] = get_ASAPfile(ID,day,sess);
                [path,videoname] = fileparts(asapfile);
                Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
                modfile = fullfile(path,[videoname,'_Mcells.mat']);
                load(modfile,'bins','onbins');
                load(Datafile,'V','F');
                disp(asapfile);
                
                Ronbins = onbins;%(bins >=1 & bins <= 8);       
                Aonbins = (V.time >=1 & V.time <= 7);
    
                f = 2;
                
                R = V.R;
                R = squeeze(R(roi,:,:));
                R = R(Ronbins,:);
                R1 = mean(R,'all');
                
                A = F.Ampl;
                A = squeeze(A(roi,:,:,f));                                   % Keep roi rates
                A = A(Aonbins,:);
                A1 = mean(A,'all');
        
                L1 = F.lVec(roi,f);
                
                for d2 = d1+1 : length(Days{i})                             % For that day and the next on the session list
                    % FIND CORRECT FILE TO LOAD
                    D2 = Days{i}(d2);
                    Sess2 = reg(D2,1);                                      % Actual imaging session number
                    day = asap(a).sessions{Sess2,1};                        % Actual date
                    sess = reg(D2,2);                                       % Actual video index
                    roi = reg(D2,3);                                        % ROI index
                    
                    % LOAD CORRECT FILE (CANNOT REPLACE WITH POOLED DATA DUE TO SESSION NUMBERING IN DAYS_PV/SOM)
                    [asapfile,~] = get_ASAPfile(ID,day,sess);
                    [path,videoname] = fileparts(asapfile);
                    Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
                    modfile = fullfile(path,[videoname,'_Mcells.mat']);
                    load(modfile,'onbins');
                    load(Datafile,'V','F');
                    disp(asapfile);
                    
                    R = V.R;
                    R = squeeze(R(roi,:,:));
                    R = R(Ronbins,:);
                    R2 = mean(R,'all');
                
                    A = F.Ampl;
                    A = squeeze(A(roi,:,:,f));                               % Keep roi rates 
                    A = A(Aonbins,:);
                    A2 = mean(A,'all');
          
                    L2 = F.lVec(roi,f);
                
                    MR{a,i} = cat(1,MR{a,i},[R1 R2]);
                    AMP{a,i} = cat(1,AMP{a,i},[A1 A2]);                         % Stack 
                    LV{a,i} = cat(1,LV{a,i},[L1 L2]);
                    Ddist{a,i} = cat(1,Ddist{a,i}, Sess2-Sess1);            % and the distance between imaging sessions
                end
            end
        end
    end
end

%% SET PRE-POST DAY DISTANCE TO 6
for a = 1:Na
    Ddist{a,2} = 6*ones(size(Ddist{a,2}));
end

%% POOL ALL DAY DISTANCES
R1 = cell2mat(MR(:,1));
A1 = cell2mat(AMP(:,1));
L1 = cell2mat(LV(:,1));
D1 = cell2mat(Ddist(:,1));

R2 = cell2mat(MR(:,2));
A2 = cell2mat(AMP(:,2));
L2 = cell2mat(LV(:,2));
D2 = cell2mat(Ddist(:,2));

%% COMPUTE DIFFERENCE AFTER - BEFORE
dR1 = diff(R1,[],2);
dA1 = diff(A1,[],2);
dL1 = diff(L1,[],2);

dR2 = diff(R2,[],2);
dA2 = diff(A2,[],2);
dL2 = diff(L2,[],2);

%% LEAST SQUARES ESTIMATE FOR DISTRIBUTION
% Rcoefficients = polyfit(D1, dR1, 1);                    % Equivalent to: [ones(size(D)), D]\C
% Acoefficients = polyfit(D1, dA1, 1);                    % Equivalent to: [ones(size(D)), D]\C
% Lcoefficients = polyfit(D1, dL1, 1);                    % Equivalent to: [ones(size(D)), D]\C
% Rfitted = polyval(Rcoefficients, [0; D1; max(D1)+1]);
% Afitted = polyval(Acoefficients, [0; D1; max(D1)+1]);
% Lfitted = polyval(Lcoefficients, [0; D1; max(D1)+1]);
% Rpfit = coefTest(fitlm(D1,dR1,'linear'));
% Apfit = coefTest(fitlm(D1,dA1,'linear'));
% Lpfit = coefTest(fitlm(D1,dL1,'linear'));

%% SPEARMAN CORRELATION OF DISTRIBUTIONS ACROSS DAY-DISTANCE
[Rr,Pr] = corr(D1,dR1,'type','Spearman');
[Ra,Pa] = corr(D1,dA1,'type','Spearman');
[Rl,Pl] = corr(D1,dL1,'type','Spearman');


%% PLOT RATES
figure; 
subplot(3,3,1:2);hold on;
plot(D1 + randn(size(D1))*0.1,dR1,'ok','MarkerFaceColor','k')
% plot([0; D1; max(D1)+1] , Rfitted,'-k')

plot(D2 + randn(size(D2))*0.1,dR2,'ok')

xlabel('Distance between sessions');
ylabel('Diff of Rates');
% title(num2str(Rpfit));
title([num2str(Rr),' (',num2str(Pr),')']);

%% PLOT AMPLITUDE 
subplot(3,3,4:5);hold on;
plot(D1 + randn(size(D1))*0.1,dA1,'ok','MarkerFaceColor','k')
% plot([0; D1; max(D1)+1] , Afitted,'-k')

plot(D2 + randn(size(D2))*0.1,dA2,'ok')

xlabel('Distance between sessions');
ylabel('Diff of Amplitude');
% title(num2str(Apfit));
title([num2str(Ra),' (',num2str(Pa),')']);

%% PLOT LVEC
subplot(3,3,7:8);hold on;
plot(D1 + randn(size(D1))*0.1,dL1,'ok','MarkerFaceColor','k')
% plot([0; D1; max(D1)+1] , Lfitted,'-k')

plot(D2 + randn(size(D2))*0.1,dL2,'ok')

xlabel('Distance between sessions');
ylabel('Diff of Vector length');
% title(num2str(Lpfit));
title([num2str(Rl),' (',num2str(Pl),')']);

%% COMPARE TRAINED VS PRE-POST
subplot(333);hold on;
p = ranksum(dR1,dR2);

cc = table(dR1',dR2','VariableNames',{'trained','pre-post'});
violinplot(cc);
plot_significance(p,1,2,dR1,dR2)
title(num2str(p));
ylabel('Rates');

%% COMPARE TRAINED VS PRE-POST
subplot(336);hold on;
p = ranksum(dA1,dA2);

cc = table(dA1',dA2','VariableNames',{'trained','pre-post'});
violinplot(cc);
plot_significance(p,1,2,dA1,dA2)
title(num2str(p));
ylabel('Amplitude');

%% COMPARE TRAINED VS PRE-POST
subplot(339);hold on;
p = ranksum(dL1,dL2);

cc = table(dL1',dL2','VariableNames',{'trained','pre-post'});
violinplot(cc);
plot_significance(p,1,2,dL1,dL2)
title(num2str(p));
ylabel('Lvec');

%% KEEP ONLY DAYS 1 SESSION APART AND PLOT  PRE-POST AND DAY X - DAY X+1
R1 = R1(D1 == 1,:); 
A1 = A1(D1 == 1,:); 
L1 = L1(D1 == 1,:); 

% RATES
figure;
subplot(311); hold on;
for i = 1:size(R2,1)
    plot(1:2,[R2(i,1), R2(i,2)],'o-','Color',[.8 .8 .8]);                    
end
[~,PR_prepost] = ttest(R2(:,1),R2(:,2))                                            
errorbar([0.9 2.1],[mean(R2(:,1)) mean(R2(:,2))],[SEM(R2(:,1)) SEM(R2(:,2))],'ok-')   
plot_significance(PR_prepost,0.9,2.1,mean(R2(:,1)),mean(R2(:,2)))

for i = 1:size(R1,1)
    plot(4:5,[R1(i,1), R1(i,2)],'o-','Color',[.8 .8 .8]);                  
end
[~,PR_post] = ttest(R1(:,1),R1(:,2))                                            
errorbar([3.9 5.1],[mean(R1(:,1)) mean(R1(:,2))],[SEM(R1(:,1)) SEM(R1(:,2))],'ok-')   
plot_significance(PR_post,3.9,5.1,mean(R1(:,1)),mean(R1(:,2)))
xlim([0.8 5.2]);

% AMPLITUDE
subplot(312); hold on;
for i = 1:size(A2,1)
    plot(1:2,[A2(i,1), A2(i,2)],'o-','Color',[.8 .8 .8]);                    
end
[~,PA_prepost] = ttest(A2(:,1),A2(:,2))                                            
errorbar([0.9 2.1],[mean(A2(:,1)) mean(A2(:,2))],[SEM(A2(:,1)) SEM(A2(:,2))],'ok-')   
plot_significance(PA_prepost,0.9,2.1,mean(A2(:,1)),mean(A2(:,2)))

for i = 1:size(A1,1)
    plot(4:5,[A1(i,1), A1(i,2)],'o-','Color',[.8 .8 .8]);                   
end
[~,PA_post] = ttest(A1(:,1),A1(:,2))                                           
errorbar([3.9 5.1],[mean(A1(:,1)) mean(A1(:,2))],[SEM(A1(:,1)) SEM(A1(:,2))],'ok-')   
plot_significance(PA_post,3.9,5.1,mean(A1(:,1)),mean(A1(:,2)))
xlim([0.8 5.2]);

% LVECTOR
subplot(313); hold on;
for i = 1:size(L2,1)
    plot(1:2,[L2(i,1), L2(i,2)],'o-','Color',[.8 .8 .8]);                    
end
[~,PL_prepost] = ttest(L2(:,1),L2(:,2))                                            
errorbar([0.9 2.1],[mean(L2(:,1)) mean(L2(:,2))],[SEM(L2(:,1)) SEM(L2(:,2))],'ok-')  
plot_significance(PL_prepost,0.9,2.1,mean(L2(:,1)),mean(L2(:,2)))

for i = 1:size(L1,1)
    plot(4:5,[L1(i,1), L1(i,2)],'o-','Color',[.8 .8 .8]);                    
end
[~,PL_post] = ttest(L1(:,1),L1(:,2))                                            
errorbar([3.9 5.1],[mean(L1(:,1)) mean(L1(:,2))],[SEM(L1(:,1)) SEM(L1(:,2))],'ok-')   
plot_significance(PL_post,3.9,5.1,mean(L1(:,1)),mean(L1(:,2)))
xlim([0.8 5.2]);

