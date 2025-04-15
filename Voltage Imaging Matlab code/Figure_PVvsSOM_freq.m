clear

Days_PV;
ASAP{1} = asap;

Days_SOM;
ASAP{2} = asap;

Fs = 1000;
window = 60;
peritime = -window/Fs : 1/Fs : window/Fs;

phases = -pi:0.2:pi;
% doublephases = [phases(1:end-1), phases(1:end-1)+2*pi];

freqtit = {'Delta (0.5-3Hz)';'Theta (4-10Hz)';'Beta (15-30Hz)';'Gamma (40-90Hz)'};

lf = 4;

%% INITIALIZE
AMP = cell(2,1);
mDir = cell(2,1);
lVec = cell(2,1);
pvalR = cell(2,1);
pvalV = cell(2,1);
K = cell(2,1);

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
                load(datafile,'F');
                
                lf = size(F.Ampl,4);
                
                for r = 1:size(F.mDir,1)
                    mDir{ct}(count,:) = F.mDir(r,:);  % Mean vector direction across frequencies
                    lVec{ct}(count,:) = F.lVec(r,:);  % Mean vector length across frequencies
                    K{ct}(count,:) = F.Skew(r,:);
                    
                    % COMPUTE BANDPOWER
                    for f = 1:lf
                        ampl = F.Ampl(r,:,:,f);
                        AMP{ct}(count,f) = mean(ampl(:));
                    end

                    count = count + 1;
                end               
            end
        end
    end
end

%% COMPARE BANDPOWER
figure;
for f = 1:lf
    subplot(4,6,(f-1)*6+5); hold on
    amp = table(AMP{1}(:,f)',AMP{2}(:,f)','VariableNames',{'PV','SOM'});
    violinplot(amp);
    pamp = ranksum(AMP{1}(:,f),AMP{2}(:,f));
    plot_significance(pamp,1,2,AMP{1}(:,f)',AMP{2}(:,f)')
    title(num2str(pamp));
    ylabel(freqtit{f});
end

%% COMPARE STRENGTH OF PHASE LOCKING VERSUS FREQUENCIES
subplot(466); hold on
for r = 1:size(lVec{1},1)
    plot(1:4,lVec{1}(r,:),'Color',[.5 .5 .5])
end
for r = 1:size(lVec{2},1)
    plot(1:4,lVec{2}(r,:),'Color',[.8 .8 .8])
end
errorbar(1:4,mean(lVec{1},1),SEM(lVec{1}),'Color','b','LineWidth',2)
errorbar(1:4,mean(lVec{2},1),SEM(lVec{2}),'Color','r','LineWidth',2)
set(gca,'Xtick',1:4,'XTickLabel',{'Delta','Theta','Beta','Gamma'});
ylabel('Mean vector length');

%% COMPARE PHASE LOCKING ANGLE VERSUS FREQUENCIES
subplot(4,6,12); hold on
for r = 1:size(mDir{1},1)
    plot(1:4,mDir{1}(r,:),'Color',[.5 .5 .5])
end
for r = 1:size(mDir{2},1)
    plot(1:4,mDir{2}(r,:),'Color',[.8 .8 .8])
end

for f = 1:4
    M1(f) = circ_mean(mDir{1}(:,f));    
    S1(f) = circ_std(mDir{1}(:,f))/sqrt(size(mDir{1},1));
    M2(f) = circ_mean(mDir{2}(:,f));    
    S2(f) = circ_std(mDir{2}(:,f))/sqrt(size(mDir{2},1));
end
errorbar(1:4,M1,S1,'Color','b','LineWidth',2)
errorbar(1:4,M2,S2,'Color','r','LineWidth',2)
set(gca,'Xtick',1:4,'XTickLabel',{'Delta','Theta','Beta','Gamma'});
ylabel('Mean phase');

%% COMPARE SKEWNESS VERSUS FREQUENCIES
subplot(4,6,18); hold on
for r = 1:size(lVec{1},1)
    plot(1:4,K{1}(r,:),'Color',[.5 .5 .5])
end
for r = 1:size(lVec{2},1)
    plot(1:4,K{2}(r,:),'Color',[.8 .8 .8])
end
errorbar(1:4,mean(K{1},1),SEM(K{1}),'Color','b','LineWidth',2)
errorbar(1:4,mean(K{2},1),SEM(K{2}),'Color','r','LineWidth',2)
set(gca,'Xtick',1:4,'XTickLabel',{'Delta','Theta','Beta','Gamma'});
ylabel('Phase skewness');


%% COMPARE PHASE LOCKING
for f = 1:lf
    subindex = (6*(4*(f-1)+1)+[1 2 3])+[0;6;12];
    subplot(4*lf,6, subindex(:));hold on;
    
    mPV = mDir{1}(:,f);                                    % Keep mean phase locking of PV
    mSOM = mDir{2}(:,f);                                    % and SOM
    lPV = lVec{1}(:,f);
    lSOM = lVec{2}(:,f);

    % Plot distributions
    plot(mPV,lPV,'vb','Markersize',3);
    plot(mSOM,lSOM,'vr','Markersize',3);
    
    % Plot oscillation signal
    mn = (min([lPV; lSOM]));                                    
    mx = (max([lPV; lSOM]));                                
    plot(phases, (cos(phases)+1)/2,'k','Linewidth',1.5);
    set(gca,'Xtick',[-pi, 0, pi, 2*pi 3*pi],'XTicklabel',{'-\pi', '0', '\pi', '2\pi', '3\pi'})
    
    if f == lf, xlabel('Phases');end
    axis tight
    box on;
    
    % Plot means and SEs of distributions
    MPV = circ_mean(mPV);
    SPV = circ_std(mPV) / sqrt(length(mPV));
    errorbar(MPV,mean(lPV),SEM(lPV),SEM(lPV),SPV,SPV,...
        's','Color','b','MarkerSize',10,'MarkerEdgeColor','b','Linewidth',2);
    
    MSOM = circ_mean(mSOM);
    SSOM = circ_std(mSOM) / sqrt(length(mSOM));
    errorbar(MSOM,mean(lSOM),SEM(lSOM),SEM(lSOM),SSOM,SSOM,...
        's','Color','r','MarkerSize',10,'MarkerEdgeColor','r','Linewidth',2);

    % Plot top histogram
    subplot(4*lf,6, 4*6*(f-1)+[1 2 3]);hold on;
    h1 = histogram(mPV,phases(1):0.05:phases(end));
    h2 = histogram(mSOM,phases(1):0.05:phases(end));
    Z = max([h1.Values, h2.Values]);
    line(MPV*[1 1],[0 Z], 'Color','b','Linewidth',2);
    line(MSOM*[1 1],[0 Z], 'Color','r','Linewidth',2);
    
    p = circ_wwtest(mPV,mSOM);                          % Compare their mean phases with parametric Watson-Williams test .
    plot_significance(p,MPV,MSOM,h1.Values,h2.Values)
    axis tight;
    title(freqtit{f});
    set(gca,'Xtick',[]);
    
    % Plot side histogram
    subplot(4*lf,6, 4*6*(f-1)+6*[2 3 4]-2);hold on; view(90, -90)
    h1 = histogram(lPV,0:0.05:0.9);
    h2 = histogram(lSOM,0:0.05:0.9);
    Z = max([h1.Values, h2.Values]);
    line(mean(lPV)*[1 1],[0 Z], 'Color','b','Linewidth',2);
    line(mean(lSOM)*[1 1],[0 Z], 'Color','r','Linewidth',2);
    
    p = ranksum(lPV,lSOM);                              % Compare their vector length medians.
    plot_significance(p,mean(lPV),mean(lSOM),h1.Values,h2.Values)
    set(gca,'Xtick',[]);
    axis tight;
end

