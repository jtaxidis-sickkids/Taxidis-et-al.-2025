clear

Days_PV;
ASAP{1} = asap;

Days_SOM;
ASAP{2} = asap;

Fs = 1000;
window = 60;
peritime = -window/Fs : 1/Fs : window/Fs;

C = [0.5 0.5 0.5; 
     0.8 0.8 0.8];

%% INITIALIZE
W = cell(2,1);
Trough = cell(2,1);
PTdist = cell(2,1);
Halfwidth = cell(2,1);
Z = cell(2,1);
Mperi = cell(2,1);
ROIarea = cell(2,1);

mR = cell(2,1);
InterSpike = cell(2,1);
BI = cell(2,1);
CSI = cell(2,1);

%% LOAD AND ORGANIZE DATA
for ct = 1:2
    asap = ASAP{ct};
    
    count = 1;

    for a = 1:length(asap)
        ID = asap(a).name;
        days = asap(a).sessions(:,1);
        ses = asap(a).sessions(:,3);
        for d = 1:length(days)
            for s = 1:length(ses{d})
                asapfile = get_ASAPfile(ID,days{d},ses{d}(s));
                disp(asapfile);
                [path,videoname] = fileparts(asapfile);
                datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
                load(datafile,'S','V');
                                
                for r = 1:length(S.Trough)
                    Trough{ct}(count) = S.Trough(r);
                    PTdist{ct}(count) = S.PTdist(r);
                    Halfwidth{ct}(count) = S.Halfwidth(r);
                    Z{ct}(count,:) = S.Points(r,:);
                    Mperi{ct}(count,:) = mean(S.PeriSp{r},1);
                    
                    mR{ct}(count) = S.mR(r);
                    BI{ct}(count) = S.BI(r);
                    CSI{ct}(count) = S.CSI(r);
                   
                    ROI = squeeze(V.ROI(r,:,:));
                    ROIarea{ct}(count) = sum(ROI(:));%bwarea(ROI)
  
                    count = count + 1;
                end
            end
        end
    end
end

%% REMOVE WAVEFORM BASELINE
NoBaselineMperi = Mperi;
for ct = 1:2
    for r = 1:size(Mperi{ct},1)
        NoBaselineMperi{ct}(r,:) = Mperi{ct}(r,:) - mean(Mperi{ct}(r,1:40));
    end
end

%% COMPUTE SPIKE AMPLITUDES
Peak = cell(1,2);
for ct = 1:2
    count = 1;
    for r = 1:size(NoBaselineMperi{ct},1)
        Peak{ct}(count) = NoBaselineMperi{ct}(r,window+2); % Peak is 1 point after middle
        count = count + 1;
    end
end

%% COMPARE WAVEFORMS
figure;
subplot(3,5,[4 5 9 10]); hold on;
for ct = 1:2
    for r = 1:size(Mperi{ct},1)
        plot(peritime,NoBaselineMperi{ct}(r,:),'Color',C(ct,:));
    end
end
fill_plot(peritime,mean(NoBaselineMperi{1},1),SEM(NoBaselineMperi{1}),'b');
fill_plot(peritime,mean(NoBaselineMperi{2},1),SEM(NoBaselineMperi{2}),'r');


%% COMPARE WAVEFORMS, ROI AREA, MEAN FIRING RATE, BURST INDEX, CSI
for ct = 1:2
    subplot(351); hold on;
    scatter3(Trough{ct},PTdist{ct},Halfwidth{ct},'filled');%,'MarkerEdgeColor',col(ct,:),'MarkerFaceColor',col(ct,:))
    xlabel('Trough');ylabel('PTdist');zlabel('Halfwidth');
    view(40,35); grid on
    
    subplot(356); hold on;
    scatter3(mR{ct},BI{ct},CSI{ct},'filled');%,'MarkerEdgeColor',col(ct,:),'MarkerFaceColor',col(ct,:))
    xlabel('Rate');ylabel('Burst Index');zlabel('CS Index');
    view(40,35); grid on
    
    subplot(3,5,11); hold on
    scatter3(ROIarea{ct},PTdist{ct},BI{ct},'filled');%,'MarkerEdgeColor',col(ct,:),'MarkerFaceColor',col(ct,:))
    xlabel('ROI area');ylabel('PTdist');zlabel('Burst Index');
    view(40,35); grid on
end

cols = cell_colors;
cols = cols([3,6],:);

subplot(352); hold on
mr = table(mR{1},mR{2},'VariableNames',{'PV','SOM'});
violinplot(mr);
pmr = ranksum(mR{1},mR{2});
plot_significance(pmr,1,2,mR{1},mR{2})
title(num2str(pmr));
ylabel('Mean rates');

subplot(357); hold on
bi = table(BI{1},BI{2},'VariableNames',{'PV','SOM'});
violinplot(bi);
pbi = ranksum(BI{1},BI{2});
plot_significance(pbi,1,2,BI{1},BI{2})
title(num2str(pbi));
ylabel('Burst index (10ms)');

subplot(3,5,12); hold on
csi = table(CSI{1},CSI{2},'VariableNames',{'PV','SOM'});
violinplot(csi);
pcsi = ranksum(CSI{1},CSI{2});
plot_significance(pcsi,1,2,CSI{1},CSI{2})
title(num2str(pcsi));
ylabel('CS Index');

%% COMPARE HALFWIDTH AND AMPLITUDE
subplot(3,5,14); hold on
hw = table(Halfwidth{1},Halfwidth{2},'VariableNames',{'PV','SOM'});
violinplot(hw);
phw = ranksum(Halfwidth{1},Halfwidth{2});
plot_significance(phw,1,2,Halfwidth{1},Halfwidth{2})
title(num2str(phw));
ylabel('Half Width');

subplot(3,5,15); hold on
pk = table(Peak{1},Peak{2},'VariableNames',{'PV','SOM'});
violinplot(pk);
ppk = ranksum(Peak{1},Peak{2});
plot_significance(ppk,1,2,Peak{1},Peak{2})
title(num2str(ppk));
ylabel('Amplitude');

