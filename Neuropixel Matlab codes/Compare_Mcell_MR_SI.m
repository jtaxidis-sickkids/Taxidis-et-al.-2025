function [SI1,SI2] = Compare_Mcell_MR_SI(MR1,MR2,F,bins,mcellflag)

%% SHAPE RATES IN PREFERRED VS NON_PREFERRED ODORS
ct = 1;

M1pref = [MR1{ct,1,1}; MR1{ct,2,2}];
M1npref = [MR1{ct,1,2}; MR1{ct,2,1}];
M2pref = [MR2{ct,1,1}; MR2{ct,2,2}]; 
M2npref = [MR2{ct,1,2}; MR2{ct,2,1}]; 

F = [F{1,1}; F{1,2}];

% COMPUTE FOR ODOR OR FOR TIME CELLS ONLY, OR ALL MCELLS
switch mcellflag
    case 'odor'
        k = (F <= 2);
    case 'delay'
        k = (F > 2);
    case'all'
        k = true(size(F));
end
F = F(k);
M1pref(~k,:) = [];
M1npref(~k,:) = [];
M2pref(~k,:) = [];
M2npref(~k,:) = [];
% -------------------------

%% GET RATES AT FIELD TIME BIN
Nm = length(F);
M1p = nan(Nm,1);
M1np = nan(Nm,1);
M2p = nan(Nm,1);
M2np = nan(Nm,1);

for c = 1:Nm                                                                % For each mcell
   k = find(bins == F(c));                                                  % find the bin index of the field time
   M1p(c) = M1pref(c,k);                                                    % Keep only rates in that bin
   M1np(c) = M1npref(c,k);                                                  
   M2p(c) = M2pref(c,k);                                                    
   M2np(c) = M2npref(c,k);                                                  
end

%% GET SI IN NON-OPTO AND OPTO TRIALS
SI1 = (M1p - M1np)./(M1p + M1np);
SI2 = (M2p - M2np)./(M2p + M2np);

% DO NOT USE ABSOLUTE VALUE BECAUSE IT IS ALREADY ODOR-AGNOSTIC (POSITIVE
% VALUE TOWARDS PREFERRED ODOR)
% SI1 = abs(SI1);       
% SI2 = abs(SI2);

%% COMPARE MEAN RATES IN FIELD WITHvsWITHOUT OPTO
figure;
subplot(231); hold on
plot([1 2],[M1p,M2p],'o-b');
errorbar([0.9 2.1],[nanmean(M1p) nanmean(M2p)],[SEM(M1p) SEM(M2p)],'ok-')   
set(gca,'Xtick',[1,2],'XTickLabel',{'no opto';'opto'});
ylabel('Mean rate');

[~,pM] = ttest(M1p,M2p)                                    
plot_significance(pM,0.9,2.1,nanmean(M1p), nanmean(M2p))

subplot(234); hold on
MRshift = (M1p - M2p);                                                      % Change of field rate over preferred trials
plot(F,MRshift,'ok');
line([1,7],mean(MRshift)*[1,1],'LineStyle','--','Color','k');
xlabel('Field time (sec)');
ylabel('Mean Rate shift');

edges = bins(1:5:end);                                                      % Make bins
de = edges(2)-edges(1);                                                     % bin size
[Mb,Sb] = bin_x_axis(F,MRshift,edges);                                      % Bin the SI distribution
fill_plot(edges(1:end-1)+de/2, Mb, Sb, 'b');                                % Plot
    
%% COMPARE SI WITHvsWITHOUT OPTO
subplot(232); hold on
plot([1 2],[SI1,SI2],'o-b');
errorbar([0.9 2.1],[nanmean(SI1) nanmean(SI2)],[SEM(SI1) SEM(SI2)],'ok-')   
set(gca,'Xtick',[1,2],'XTickLabel',{'no opto';'opto'});
ylabel('SI');

[~,pSI] = ttest(SI1,SI2)                                    
plot_significance(pSI,0.9,2.1,nanmean(SI1), nanmean(SI2))

subplot(235); hold on
SIshift = abs(SI1 - SI2);                                                   % Absolute change of SI
plot(F,SIshift,'ok');
line([1,7],mean(SIshift)*[1,1],'LineStyle','--','Color','k');
xlabel('Field time (sec)');
ylabel('|SI shift|');

edges = bins(1:10:end);                                                      % Make bins
de = edges(2)-edges(1);                                                     % bin size
[Mb,Sb] = bin_x_axis(F,SIshift,edges);                                      % Bin the SI distribution
fill_plot(edges(1:end-1)+de/2, Mb, Sb, 'b');                                % Plot

if strcmp(mcellflag,'all')
    subplot(236); hold on
    SIshift_od = SIshift(F <= 2);
    SIshift_de = SIshift(F > 2);
    pSI_od_de = ones(2,2);
    pSI_od_de(1,2) = ranksum(SIshift_od,SIshift_de)
    plot_mean_SE({SIshift_od,SIshift_de}, pSI_od_de, rand(2,3))
    set(gca,'Xtick',[1,2],'XTickLabel',{'odor';'delay'});
    ylabel('|SI shift|');
end


