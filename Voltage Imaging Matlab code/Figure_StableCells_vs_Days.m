function [ODORSTABILITY,INFLOW,OUTFLOW] = Figure_StableCells_vs_Days(INclass)

if strcmp(INclass,'PV')
    Days_PV
    maxd = 7;
elseif strcmp(INclass,'SOM')
    Days_SOM
    maxd = 8;
end

Na = length(asap);                                                          % total number of animals
Nr = length([asap.samecells]');                                             % total number of cells followed across days

CT = nan(Nr,maxd);
MB = nan(Nr,maxd);
SI = nan(Nr,maxd);

ISI = nan(Nr,maxd);
MR = nan(Nr,maxd);
BI = nan(Nr,maxd);

AMP = nan(Nr,maxd,4);
LVC = nan(Nr,maxd,4);

count = 1;

%% PLOT RATES AND FOV OF REGISTERED CELLS ACROSS DAYS
for a = 1:Na
    ID = asap(a).name;
    ls = length(asap(a).samecells);
    
    for s = 1:ls                                                            % For each group of registered cells
        reg = asap(a).samecells{s};
        lr = size(reg,1);                                                   % # days a cell was followed
        for r = 1:lr                                                        % For each ROI (day) within this group
            % FIND CORRECT FILE TO LOAD
            dd = reg(r,1);                                                  % Index of imaging day
            day = asap(a).sessions{dd,1};                                   % Actual date
            sess = reg(r,2);                                                % Actual session index
            roi = reg(r,3);                                                 % ROI index
            
            % LOAD CORRECT FILE (CANNOT REPLACE WITH POOLED DATA DUE TO SESSION NUMBERING IN DAYS_PV/SOM)
            [asapfile,~] = get_ASAPfile(ID,day,sess);
            [path,videoname] = fileparts(asapfile);
            Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
            modfile = fullfile(path,[videoname,'_Mcells.mat']);
            load(modfile,'Mcells','onbins');
            load(Datafile,'V','F','S');
            disp(asapfile);
            
            % KEEP CELL TYPE, FIELD AND SI
            mcells = Mcells(roi,:);                                         % Load Mcell identity of roi
            CT(count,dd) = mcells(1);                                       % Keep the cell type
            MB(count,dd) = mcells(2);                                       % and its field bin
            SI(count,dd) = mcells(4);                                       % and its selectivity index
    
            % KEEP MEAN ISI
            isi = S.ISI(roi,:);
            isi = cellfun(@nanmean, isi);
            ISI(count,dd) = nanmean(isi);
            
            % KEEP MEAN AMP AND LVEC
            amp = squeeze(F.Ampl(roi,onbins,:,:));
            AMP(count,dd,:) = nanmean(amp,[1,2]);
            
            LVC(count,dd,:) = F.lVec(roi,:);
            
            % KEEP MEAN RATE
            rr = V.R(roi,onbins,:);
            MR(count,dd) = mean(rr,'all');

            % KEEP BI
            BI(count,dd) = S.BI(roi);
        end
        count = count + 1;
        
        disp(' ');
    end
end

%% PLOT EVOLUTION OF FIELDS AND SI
figure;
for r = 1:Nr                                                                % For each registered cell
    days = find(~isnan(MB(r,:)));                                           % Find all the days it was imaged
    isA = find(CT(r,:) == 1);                                               % Find days where it's odorA-selective
    isB = find(CT(r,:) == 2);                                               % Find days where it's odorB-selective
    isM = find(CT(r,:) == 0);                                               % Find days where it's non-odor-selective
    
    subplot(211); hold on; view(90, -90);set(gca,'Xdir','reverse');
    plot(days,MB(r,days),'-ob');
    plot(isA,MB(r,isA),'oy','MarkerFaceColor','y');
    plot(isB,MB(r,isB),'og','MarkerFaceColor','g');
    plot(isM,MB(r,isM),'ob','MarkerFaceColor','b');
    view(90, -90)
    
    subplot(212); hold on;
    plot(days,SI(r,days),'-ob');
    plot(isA,SI(r,isA),'oy','MarkerFaceColor','y');
    plot(isB,SI(r,isB),'og','MarkerFaceColor','g');
    plot(isM,SI(r,isM),'ob','MarkerFaceColor','b');
    ylabel('SI');
end

%% PLOT EVOLUTION OF mR, BI, ISI
figure;
for r = 1:Nr                                                                % For each registered cell
    days = find(~isnan(MB(r,:)));                                           % Find all the days it was imaged
%     isA = find(CT(r,:) == 1);                                               % Find days where it's odorA-selective
%     isB = find(CT(r,:) == 2);                                               % Find days where it's odorB-selective
%     isM = find(CT(r,:) == 0);                                               % Find days where it's non-odor-selective
    
    subplot(341); hold on;
    plot(days,ISI(r,days),'-ob');
%     plot(isA,ISI(r,isA),'oy','MarkerFaceColor','y');
%     plot(isB,ISI(r,isB),'og','MarkerFaceColor','g');
%     plot(isM,ISI(r,isM),'ob','MarkerFaceColor','b');
    ylabel('ISI');
    
    subplot(342); hold on;
    plot(days,MR(r,days),'-ob');
%     plot(isA,MR(r,isA),'oy','MarkerFaceColor','y');
%     plot(isB,MR(r,isB),'og','MarkerFaceColor','g');
%     plot(isM,MR(r,isM),'ob','MarkerFaceColor','b');
    ylabel('mean rate');
    
    subplot(343); hold on;
    plot(days,BI(r,days),'-ob');
%     plot(isA,BI(r,isA),'oy','MarkerFaceColor','y');
%     plot(isB,BI(r,isB),'og','MarkerFaceColor','g');
%     plot(isM,BI(r,isM),'ob','MarkerFaceColor','b');
    ylabel('BI');
end

%% PLOT EVOLUTION OF AMPLITUDES AND VECTOR LENGTHS
for r = 1:Nr 
    days = find(~isnan(MB(r,:)));                                           % Find all the days it was imaged
%     isA = find(CT(r,:) == 1);                                               % Find days where it's odorA-selective
%     isB = find(CT(r,:) == 2);                                               % Find days where it's odorB-selective
%     isM = find(CT(r,:) == 0);                                               % Find days where it's non-odor-selective

    for f = 1:4
        subplot(3,4,4+f); hold on;
        plot(days,AMP(r,days,f),'-ob');
%         plot(isA,AMP(r,isA,f),'oy','MarkerFaceColor','y');
%         plot(isB,AMP(r,isB,f),'og','MarkerFaceColor','g');
%         plot(isM,AMP(r,isM,f),'ob','MarkerFaceColor','b');
        ylabel('Amplitude');
        
        subplot(3,4,8+f); hold on;
        plot(days,LVC(r,days,f),'-ob');
%         plot(isA,LVC(r,isA,f),'oy','MarkerFaceColor','y');
%         plot(isB,LVC(r,isB,f),'og','MarkerFaceColor','g');
%         plot(isM,LVC(r,isM,f),'ob','MarkerFaceColor','b');
        ylabel('Vector length');
    end
end

%% PLOT FLUX INTO ODOR-FIELDS
odorfield = double(MB <= 2);
odorfield(isnan(MB)) = nan;
FLOW = nan(Nr,maxd);
ODORSTABILITY = nan(Nr,maxd);
for r = 1:Nr
    k = find(~isnan(odorfield(r,:)));                                       % Keep indexes of days the cell was recorded
    for i = 2:length(k)                                                     % For each recorded day (excpect first)
        FLOW(r,k(i)) = odorfield(r,k(i)) - odorfield(r,k(i-1));
        ODORSTABILITY(r,k(i)) = (odorfield(r,k(i)) == 1 & odorfield(r,k(i-1)) == 1);
    end
end

INFLOW = 100 * nansum(FLOW==1,1) ./ sum(~isnan(FLOW),1);
OUTFLOW = 100 * nansum(FLOW==-1,1) ./ sum(~isnan(FLOW),1);
ODORSTABILITY = 100 * nansum(ODORSTABILITY,1) ./ sum(~isnan(ODORSTABILITY),1);

% PLOT ACROSS DAYS
figure;
subplot(1,3,1:2);
plot(1:maxd,INFLOW,1:maxd,OUTFLOW,1:maxd,ODORSTABILITY,'o-')

% PLOT AVERAGE AND COMPARE
subplot(133);hold on;
p = nan(3,1);
p(1) = ranksum(ODORSTABILITY,INFLOW);
p(2) = ranksum(ODORSTABILITY,OUTFLOW);
p(3) = ranksum(INFLOW,OUTFLOW);
[~,p] = fdr(p);

cc = table(ODORSTABILITY',INFLOW',OUTFLOW','VariableNames',{'stable','inflow','outflow'});
violinplot(cc);
plot_significance(p(1),1,2,ODORSTABILITY,INFLOW);
plot_significance(p(2),1,3,ODORSTABILITY,OUTFLOW);
plot_significance(p(3),2,3,INFLOW,OUTFLOW);
title(num2str(p));
ylabel('Mean flow across days');
