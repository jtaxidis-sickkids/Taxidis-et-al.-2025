function Figure_DFF_stats(INclass)  

%% LOAD POOLED DATA AND STACK
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'DFF','MC','TR','timepoints','time');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    DFF{1,5} = [];     % DFF{2,3} = [];
    MC{1,5} = [];      % MC{2,3} = [];
    TR{1,5} = [];      % TR{2,3} = [];
end
DFF = vertcat(DFF{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

ls = length(DFF);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
DFFall = {};
MCall = [];
TRall = {};

for i = 1:ls
    for r = 1:size(DFF{i},1)
        DFFall{count,1} = squeeze(DFF{i}(r,:,:))';                          % Keep only that cell's DFF [trials x time]
        MCall(count,:) = MC{i}(r,:);    
        TRall{count,1} = TR{i};
        count = count + 1;
    end
end
DFF = DFFall;
MC = MCall;
TR = TRall;
clear DFFall MCall TRall

ls = count-1;

%% MAX AND MIN DFF, DIP ONSET TIME AND DURATION IN EACH TRIAL
MinD = cell(ls,1);   
MinT = cell(ls,1);   
MaxD = cell(ls,1);   
MinOn = cell(ls,1);  
MinDur = cell(ls,1);  
sigDip = cell(ls,1);

for i = 1:ls                                                                % For each cell
    dff = sgolayfilt(DFF{i}, 1, 21, [] , 2);                                % Apply linear filter along time bins
    
    [MinD{i},MaxD{i},oddff] = get_DFF_minmax(dff,time,timepoints);          % NORMALIZE ODOR DFF TO BASELINE AND COMPUTE MIN/MAX

    Ntr = size(oddff,1);
    MinOn{i} = nan(Ntr,1);
    MinDur{i} = nan(Ntr,1);
    MinT{i} = nan(Ntr,1);
    
    for tr = 1:Ntr                                                          % For each trial      
        tdip_low = find(oddff(tr,:) < -2);                                  % find any zscored odor DFFs less than -2
%         figure;hold on; c = 0;
        while any(diff(tdip_low) > 1)                                       % while there are disconnected low-DFF segments                             
%             plot(tdip_low,c*ones(size(tdip_low)),'o-');c=c+1;
            contin_dip = diff(tdip_low) == 1;                               % find points corresponding to continuous segments
            tdip_low = tdip_low(contin_dip);                                % keep only those
        end                                                                 % until there is only one continuous segment
        
        if length(tdip_low) > 50                                            % If that segment is at least 50msec long
            tdip_high = find(oddff(tr,tdip_low) < -4);                      % Find if any dff sub-segments lower than 4  
            tdip_high = tdip_high + tdip_low(1) - 1;                        % keep the actual indexes of those segments 
            
            while any(diff(tdip_high) > 1)                                  % Repeat to keep only one continuous segment
                contin_dip = diff(tdip_high) == 1;
                tdip_high = tdip_high(contin_dip);
            end
            
            if length(tdip_high) > 10                                       % If that segment is at least 10msec long
                [~,tmin] = min(oddff(tr,tdip_high));                        % Keep the timing of the minimum
                MinT{i}(tr) = tmin/1000 + tdip_high(1)/1000 + timepoints(1); % turn it to actual timepoint
                
                MinDur{i}(tr) = (tdip_low(end) - tdip_low(1)) / 1000;       % Keep duration of lower-thresh segment (in sec)
                MinOn{i}(tr) = tdip_low(1)/1000 + timepoints(1);            % Keep initial point of lower-thresh segment
            end
        end
    end

    % Count "SIGNIFICANT" DIPS
    sigDip{i} = ~isnan(MinT{i});
end

%% SPLIT MINIMUMS OF CELL GROUPS
mcMin = cell(4,1);
mcTR = cell(4,1);
mcSig = cell(4,1);

MC = MC(:,1);                                                               % Keep only cell types

for tt = 1:4
    k = find(MC == tt);
    if tt == 3
        k = find(MC == 0);
    elseif tt == 4
        k = find(isnan(MC));
    end
    mcMin{tt} = MinD(k);
    mcTR{tt} = TR(k);
    mcSig{tt} = sigDip(k);
end

%% HISTOGRAM OF MINS AND MAXS ACROSS ALL TRIALS AND THEIR RELATIONSHIP
mind = cell2mat(MinD);                                                      % Pool all cells and all trials
maxd = cell2mat(MaxD);

% PLOT ONLY SIGNIFICANT MINS AND MAXS
% sig = cell2mat(sigDip);
% mind=mind(sig);
% maxd=maxd(sig);

[mean(mind) std(mind)]

figure;
subplot(331); hold on;
histogram(mind,[-10:0.05:0]);
line(nanmedian(mind)*[1 1],[0 50]);
xlabel('Min DFF'); ylabel('#trials');

subplot(332); hold on;
histogram(maxd,[0:0.05:10])
line(nanmedian(maxd)*[1 1],[0 45]);
xlabel('Max DFF'); ylabel('#trials');

% PLOT RELATIONSHIP BETWEEN MIN AND MAX PER TRIAL
subplot(333); hold on
scatter(abs(mind),maxd);
xlabel('|Min DFF|'); ylabel('Max DFF');

% LEAST SQUARES ESTIMATE FOR DISTRIBUTION
coefficients = polyfit(abs(mind),maxd, 1);                    % Equivalent to: [ones(size(D)), D]\C
Cfitted = polyval(coefficients, [min(abs(mind)) max(abs(mind))]);
pfit = coefTest(fitlm(abs(mind),maxd,'linear'));
plot([min(abs(mind)) max(abs(mind))] , Cfitted,'-k')
title(num2str(pfit));

%% DISTRIBUTION OF AVERAGE MINS PER CELL
mind = cellfun(@mean, MinD);                                                % Average minimum per cell across all trials

subplot(334); hold on;
histogram(mind,[-6:0.05:-2])
line(nanmedian(mind)*[1 1],[0 7]);
xlabel('<Min DFF>'); ylabel('#cells');

%% PLOT DISTRIBUTION OF % SIGNIFICANT DIP TRIALS PER CELL
dip = cellfun(@sum,sigDip) ./  cellfun(@length,sigDip) * 100;                     % Ratio of significant dips per cell

subplot(335); hold on;
histogram(dip,[0:2:100])
line(nanmedian(dip)*[1 1],[0 30]);
xlabel('% signif dips'); ylabel('#cells');

%% PLOT DISTRIBUTION OF TIMING OF DIPS
tmin = cell2mat(MinT);                                              % Timing of minimum (in sec)

subplot(336); hold on;
histogram(tmin,[1:0.01:2])
line(nanmedian(tmin)*[1 1],[0 20]);
xlabel('timing of Min'); ylabel('#trials');

%% PLOT DISTRIBUTION OF DURATION OF DIPS
mindur = cell2mat(MinDur);                                                      % Pool all cells and all trials

subplot(337); hold on;
histogram(mindur,[0:0.01:1])
line(nanmedian(mindur)*[1 1],[0 20]);
xlabel('dip duration'); ylabel('#trials');

%% PLOT DISTRIBUTION OF ONSET OF DIPS
minon = cell2mat(MinOn);                                                      % Pool all cells and all trials

subplot(338); hold on;
histogram(minon,[1:0.01:2])
line(nanmedian(minon)*[1 1],[0 20]);
xlabel('dip onset (sec)'); ylabel('#trials');

%%
nanmedian(minon)
nanmedian(tmin)
nanmedian(mindur)

%% PLOT <DFF> ABOVE VERSUS BELOW CRITERION
D1 = [];
D2 = [];
for i = 1:ls
    dip = sigDip{i};                                                           % Keep the significant-dip indexes 
    D1 = cat(1, D1, DFF{i}(dip,:));                                         % Keep corresponding DFFs
    D2 = cat(1, D2, DFF{i}(~dip,:));
end

subplot(339); hold on;
fill_plot(time,nanmean(D1,1),SEM(D1),'k');
fill_plot(time,nanmean(D2,1),SEM(D2),'r');


%% PLOT DIP IN PREFERRED - NONPREFERRED IN ODOR-SPECIFIC MCELLS
PREF = [];
NPREF = [];
for tt = 1:2                                                                % For each odot
   for i = 1:sum(MC == tt)                                                  % For each cell specific to that odor
        mins = mcMin{tt}{i};                                                % Keep is minima over all trials
        trs = mcTR{tt}{i};                                                  % and the corresponding trials
        
        PREF = [PREF; mean(mins(trs == tt))];                               % Store average minimum of preferred trials
        NPREF = [NPREF; mean(mins(trs~=tt))];                               % and non preferred
   end
end
p = ranksum(PREF,NPREF);                                                    % Compare

figure;
subplot(221); hold on
for i = 1:length(PREF)
    plot(1:2,[PREF(i) NPREF(i)],'o-','Color',[0.8 0.8 0.8]);                    % Plot mean ON vs mean OFF rates
end
errorbar([0.9 2.1],[mean(PREF) mean(NPREF)],[SEM(PREF) SEM(NPREF)],'ok-')   % Plot means over all cells
plot_significance(p,0.9,2.1,mean(PREF),mean(NPREF))
xlim([0.8 2.2]);
ylabel('Pooled Min');
set(gca,'XTick',[1 2],'XTickLabel',{'Preferred','Non-preferred'});

%% COMPARE DIP BETWEEN GROUPS
mMin = cell(4,1);
for i = 1:4
    mMin{i} = cellfun(@mean, mcMin{i});
end

p = ones(1,3);
p(1) = ranksum([mMin{1}; mMin{2}],mMin{3});                                 % Compare ALL odor-specific with non-odor-specific
p(2) = ranksum([mMin{1}; mMin{2}],mMin{4});                                 % and with nonMcells
p(3) = ranksum(mMin{3},mMin{4});                                        

subplot(222); hold on
hw = table([mMin{1}; mMin{2}]',mMin{3}',mMin{4}','VariableNames',{'Odor-Sp','Non-Odor-Sp','No Field'});
violinplot(hw);
plot_significance(p(1),1,2,[mMin{1}; mMin{2}],mMin{3})
plot_significance(p(2),1,2,[mMin{1}; mMin{2}],mMin{4})
plot_significance(p(3),1,2, mMin{3},mMin{4})
title(num2str(p));
ylabel('<Min> over trials');

%% COMPARE # SIGNIFICANT DIPS BETWEEN GROUPS
msig = cell(4,1);
for i = 1:4
    msig{i} = cellfun(@sum, mcSig{i}) ./ cellfun(@length, mcSig{i}) * 100;
end

p = ones(1,3);
p(1) = ranksum([msig{1}; msig{2}],msig{3});                                 % Compare ALL odor-specific with non-odor-specific
p(2) = ranksum([msig{1}; msig{2}],msig{4});                                 % and with nonMcells
p(3) = ranksum(msig{3},msig{4});                                        

subplot(224); hold on
hw = table([msig{1}; msig{2}]',msig{3}',msig{4}','VariableNames',{'Odor-Sp','Non-Odor-Sp','No Field'});
violinplot(hw);
plot_significance(p(1),1,2,[msig{1}; msig{2}],msig{3})
plot_significance(p(2),1,2,[msig{1}; msig{2}],msig{4})
plot_significance(p(3),1,2, msig{3},msig{4})
title(num2str(p));
ylabel('% Sig Dips');
