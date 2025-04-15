clear

Days_PV;
ASAP{1} = asap;

Days_SOM;
ASAP{2} = asap;

%% INITIALIZE
NR = cell(2,1);
SI = cell(2,1);
SI_RAND = cell(2,1);
LMC = cell(2,1);
MC = cell(2,1);
OD = cell(2,1); 
DE = cell(2,1); 

%% LOAD AND ORGANIZE DATA
for ct = 1:2
    asap = ASAP{ct};

    si = [];
    SI_rand = [];
    Lmc = [];
    Nr = [];
    MC{ct} = {};
    Od = [];
    De = [];
    
    for a = 1:length(asap)                                                  % For each animal
        ID = asap(a).name;
        days = asap(a).sessions(:,1);
        ses = asap(a).sessions(:,3);
        for d = 1:length(days)        % CHANGE TO 1:2 FOR NAIVE ONLY OR TO 3:LENGTH(DAYS) FOR TRAINED ONLY (REVISIONS)
            lmc = 0;
            nr = 0;
            od = 0;
            de = 0;
            
            for s = 1:length(ses{d})                                        % For each video
                asapfile = get_ASAPfile(ID,days{d},ses{d}(s));
                disp(asapfile);
                [path,videoname] = fileparts(asapfile);
                modfile = fullfile(path,[videoname,'_Mcells.mat']);
                load(modfile,'Mcells','trials','bins','timepoints');
                
                MC{ct} = cat(1,MC{ct},Mcells);                              % Concatenate all cells
                
                ismcells = ~isnan(Mcells(:,1));                             % Index of Mcells
                lmc = lmc + sum(ismcells);                                  % Count Mcells
                nr = nr + length(Mcells(:,1));                              % Count all cells
                od = od + sum(Mcells(ismcells,2) <= timepoints(2));         % 
                de = de + sum(Mcells(ismcells,2) >  timepoints(2));         % 
                si = [si, Mcells(:,4)'];                                    % Concatenate SI of ALL cells
                
                % MAKE RANDOM SI DISTRIBUTION (REVISIONS)
                modfile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
                load(modfile,'V');
     
                for c = 1:length(Mcells(:,1))                               % For each cell in the video
                    f = find(bins == Mcells(c,2));                          % Find its activation peak time
                    r = squeeze(V.R(c,f,:));                                % Keep rates there for all trials
                    Ntr = length(r);                        
                    
                    si_rand = nan(1,1000);                                  % Make a vector for randomized SI
                    for rep = 1:1000                                        % For each repetition
                        k = randi([0, 1], [Ntr, 1]);                        % Make a vector of randomly picked trials
                        k = logical(k);
                        r1 = mean(r(k));                                    % Mean rate in these trials
                        r2 = mean(r(~k));                                   % and in the remaining trials
                        si_rand(rep) = (r1 - r2) / (r1 + r2);               % randomized SI
                    end
                    SI_rand = [SI_rand; si_rand];
                end
                % --------------------------------------
                
            end
            Lmc = [Lmc, lmc];
            Nr = [Nr, nr];
            Od = [Od, od];
            De = [De, de];
        end
    end
    
    LMC{ct} = Lmc;
    SI{ct} = abs(si);
    SI_RAND{ct} = abs(SI_rand);   
    NR{ct} = Nr;                
    OD{ct} = Od;
    DE{ct} = De;
    
    SI{ct}(SI{ct} > 5) = [];                                                % IT REMOVES ONE CRAZY OUTLIER!!!!!!
end

%% COUNT PERCENTAGES OF MCELLS, ODOR-CELLS, DELAY-CELLS
for ct = 1:2
    OD{ct} = OD{ct} ./ LMC{ct} * 100;
    DE{ct} = DE{ct} ./ LMC{ct} * 100;
    
    LMC{ct} = LMC{ct}./NR{ct}*100;
end

%% COUNT PERCENTAGES OF CELL CLASSES
od1MC = nan(2,1);
od2MC = nan(2,1);
nodMC = nan(2,1);
nMC = nan(2,1);
Nr = nan(2,1);

for ct = 1:2
    mc = cell2mat(MC{ct});
    mc = mc(:,1);
    od1MC(ct) = sum(mc == 1);
    od2MC(ct) = sum(mc == 2);
    nodMC(ct) = sum(mc == 0);
    nMC(ct) = sum(isnan(mc));
    Nr(ct) = length(mc);
end
od1MC = od1MC ./ Nr * 100;
od2MC = od2MC ./ Nr * 100;
nodMC = nodMC ./ Nr * 100;
nMC = nMC ./ Nr * 100;

%% COMPARE %ODORvsDELAY CELLS
p = nan(1,2);
for ct = 1:2
    [~,p(ct)] = ttest(OD{ct},DE{ct});                                       % Paired sample ttest
end

figure
for ct = 1:2
    subplot(2,1,ct); hold on;
    
    scatter_plot_comparison(DE{ct} + randn(size(DE{ct})), OD{ct} + randn(size(OD{ct})), p(ct));
    xlabel('% delay cells');
    ylabel('% odor cells');
    title(num2str(p(ct)));
end

%% COMPARE NUMBER, SI, ACTIVATION PROB AND VARIANCE OF MCELLS
figure;
subplot(221); hold on
lmc = table(LMC{1},LMC{2},'VariableNames',{'PV','SOM'});
violinplot(lmc);
plmc = ranksum(LMC{1},LMC{2});
plot_significance(plmc,1,2,LMC{1},LMC{2})
title(num2str(plmc));
ylabel('%Mcells');

subplot(223); hold on
nr = table(NR{1},NR{2},'VariableNames',{'PV','SOM'});
violinplot(nr);
pnr = ranksum(NR{1},NR{2});
plot_significance(pnr,1,2,NR{1},NR{2})
title(num2str(pnr));
ylabel('Nr');

subplot(224); hold on
si = table(SI{1},SI{2},'VariableNames',{'PV','SOM'});
violinplot(si);
psi = ranksum(SI{1},SI{2});
plot_significance(psi,1,2,SI{1},SI{2})
title(num2str(psi));
ylabel('SI');

SI_RAND{1} = SI_RAND{1}(:);
SI_RAND{2} = SI_RAND{2}(:);
errorbar(0.8,nanmean(SI_RAND{1}), SEM(SI_RAND{1}),'Color','k')
errorbar(2.2,nanmean(SI_RAND{2}), SEM(SI_RAND{2}),'Color','k')

[~,p_si(1)] = ttest(SI{1},nanmean(SI_RAND{1}));
[~,p_si(2)] = ttest(SI{2},nanmean(SI_RAND{2}));
p_si

%% PLOT PERCENTAGES OF EACH CELL CLASS
subplot(222); hold on
bar([nMC, nodMC, od2MC, od1MC],'stacked')
xlim([0.4 2.6])
set(gca,'Xtick',[1 2],'XTickLabel',{'PV';'SOM'});
