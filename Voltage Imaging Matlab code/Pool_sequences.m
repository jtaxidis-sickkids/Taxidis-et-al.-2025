function [lmc, Nr] = Pool_sequences(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
elseif strcmp(INclass,'SOM')
    Days_SOM;
end

Fp = cell(1,2);
MR = cell(1,2);
MR0 = [];
nMR = [];
MRall = [];

lmc = 0;
Nr = 0;

cols = cell_colors;

%% LOAD POOLED DATA AND STACK IN SINGLE COLUMNS
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'R','MC','TR','timepoints','bins','onbins');

% USE ONLY NAIVE OR TRAINED (REVISIONS) -------------------------
% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
% if strcmp(INclass,'SOM')
%     R(3,3:7) = R(3,4:8);
%     MC(3,3:7) = MC(3,4:8);
%     TR(3,3:7) = TR(3,4:8);
% end
% 
% R(:,1:2) = [];
% MC(:,1:2) = [];
% TR(:,1:2) = [];
% --------------------

R = vertcat(R{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

ls = length(R);                                                             % Number of sessions

%% COUNT CELLS
for i = 1:ls
    lmc = lmc + sum(~isnan(MC{i}(:,1)));                               % Count mcells
%     lmc = lmc + sum(~isnan(MC{i}(:,1)) & MC{i}(:,2) > 2);% & MC{i}(:,2) <4);   % Count delay cells
    Nr = Nr + size(MC{i},1); 
end

%% ZSCORE RATES
for i = 1:ls
    R{i} = (R{i} - mean(R{i}(:,onbins,:),2)) ./ std(R{i}(:,onbins,:),[],2); % ZSCORE THE RATE OVER MODULATION BINS ONLY!!!
    R{i}(isnan(R{i}) | isinf(R{i})) = 0;
end

%% AVERAGE MCELL RATES
for i = 1:ls
    for tt1 = 1:2                                                           % For each first odor
        k = (MC{i}(:,1) == tt1);
        mr = R{i}(k,:,:);
        mr1 = mr(:,:,TR{i} == 1);
        mr2 = mr(:,:,TR{i} == 2);
        mr1 = mean(mr1,3);
        mr2 = mean(mr2,3);
        
        mr = cat(3,mr1,mr2);
        MR{tt1} = [MR{tt1};mr];
    end
end

%% Repeat for non-odor-specific Mcells
for i = 1:ls
    k = (MC{i}(:,1) == 0);
    mr = R{i}(k,:,:);
    mr1 = mr(:,:,TR{i} == 1);
    mr2 = mr(:,:,TR{i} == 2);
    mr1 = mean(mr1,3);
    mr2 = mean(mr2,3);
    mr3 = mean(mr,3);
    
    mr = cat(3,mr1,mr2,mr3);
    MR0 = [MR0;mr];
end

%% Repeat for ALL Mcells
for i = 1:ls
    k = ~isnan(MC{i}(:,1));
    mr = R{i}(k,:,:);
    mr1 = mr(:,:,TR{i} == 1);
    mr2 = mr(:,:,TR{i} == 2);
    mr1 = mean(mr1,3);
    mr2 = mean(mr2,3);
    mr3 = mean(mr,3);
    
    mr = cat(3,mr1,mr2,mr3);
    MRall = [MRall;mr];
end

%% REPEAT FOR NON-MCELLS
for i = 1:ls 
    k = isnan(MC{i}(:,1));
    mr = R{i}(k,:,:);
    mr1 = mr(:,:,TR{i} == 1);
    mr2 = mr(:,:,TR{i} == 2);
    mr1 = mean(mr1,3);
    mr2 = mean(mr2,3);
    mr3 = mean(mr,3);
    
    mr = cat(3,mr1,mr2,mr3);
    nMR = [nMR;mr];
end
    
%% CONCATENATE AND SORT 
MC = cell2mat(MC);

for tt1 = 1:2                                                               % For each first odor
    k = (MC(:,1) == tt1);    
    F = MC(k,2);                                                % Concatenate fields
    [F,order] = sort(F);                                            % Sort fields
    Fp{tt1} = F;                                               % and fields
    MR{tt1} = MR{tt1}(order,:,:);
end

k = (MC(:,1) == 0);
F = MC(k,2);                                                % Concatenate fields
[Fp0,order] = sort(F);                                            % Sort fields
MR0 = MR0(order,:,:);
  
k = ~isnan(MC(:,1));
F = MC(k,2);                                                % Concatenate fields
[Fall,order] = sort(F);                                            % Sort fields
MRall = MRall(order,:,:);
  
k = isnan(MC(:,1));
F = MC(k,2);                                                % Concatenate fields
[nFp,order] = sort(F);                                            % Sort fields
nMR = nMR(order,:,:);

%% PLOT AVERAGE FIRING RATES OF ODOR-SPECIFIC MCELLS
figure('Name','Pooled M-cells');
for tt1 = 1:2                                                  % For each trial type
    for tt2 = 1:2                                              % For each trial type (again)
        subplot(4,4,(tt1-1)*4 + tt2); hold on;
        plot_seq_rates(MR{tt1}(:,:,tt2),[],bins,timepoints,1); % Plot pooled sequences
        plot(Fp{tt1},1:length(Fp{tt1}),'.k')
        title('Trials by first odor');
        xlim([0, timepoints(5)]);
    end
end


%% PLOT ALL ODOR-SPECIFIC MCELLS ON THEIR PREFERRED VS NON-PREFERRED TRIALS
F = cell2mat(Fp');
[F,order] = sort(F);

Mp = [MR{1}(:,:,1); MR{2}(:,:,2)];
Mp = Mp(order,:);

subplot(4,4,[3 7]);
plot_seq_rates(Mp,[],bins,timepoints,1);       
title('Preferred trials');
xlim([0, timepoints(5)]);

Mnp = [MR{1}(:,:,2); MR{2}(:,:,1)];
Mnp = Mnp(order,:);

subplot(4,4,[4 8]);
plot_seq_rates(Mnp,[],bins,timepoints,1);      
title('Non-preferred trials');
xlim([0, timepoints(5)]);

%% PLOT NON-ODOR-SPECIFIC MCELLS ON ODOR1 VS ODOR2 TRIALS AND ON ALL TRIALS
for tt = 1:2
    subplot(4,4,8+[tt tt+4]);
    plot_seq_rates(MR0(:,:,tt),[],bins,timepoints,1);      
    plot(Fp0,1:length(Fp0),'.k')
    xlim([0, timepoints(5)]);
    title(['Odor ',num2str(tt)]);
end

%% PLOT ALL MCELLS AND ALL NON-MCELLS ON ODOR1 VS ODOR2 TRIALS AND ON ALL TRIALS
figure('Name','ALL M-cells');
for tt = 1:3
    subplot(8,3,(tt-1)+[1 4 7]);
    plot_seq_rates(MRall(:,:,tt),[],bins,timepoints,1);       
    plot(Fall,1:length(Fall),'.k')
    xlim([0, timepoints(5)]);
    
    subplot(8,3,(tt-1)+10); hold on;
    plot_mrate_traces(mean(MRall(:,:,tt),1),SEM(MRall(:,:,tt)),bins,timepoints,cols(tt,:));
    
    % ------
    
    subplot(8,3,(tt-1)+[1 4 7] + 12);
    plot_seq_rates(nMR(:,:,tt),[],bins,timepoints,1);       
    plot(nFp,1:length(nFp),'.k')
    xlim([0, timepoints(5)]);
    
    subplot(8,3,(tt-1)+10 + 12); hold on;
    plot_mrate_traces(mean(nMR(:,:,tt),1),SEM(nMR(:,:,tt)),bins,timepoints,cols(tt,:));
end


