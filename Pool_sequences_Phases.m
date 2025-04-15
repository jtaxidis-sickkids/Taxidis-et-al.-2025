function Pool_sequences_Phases(INclass)

addpath('CircStat2012a');

%% INITIALIZE
cols = cell_colors;
freqtit = {'Delta (0.5-3Hz)';'Theta (4-10Hz)';'Beta (15-30Hz)';'Gamma (40-90Hz)'};

Fs = 1000;
time = 0 : 1/Fs : (11*Fs-1)/Fs;

%% LOAD POOLED DATA AND STACK IN SINGLE COLUMN
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'PH','MC','TR','timepoints');

PH = vertcat(PH{:});
MC = vertcat(MC{:});
TR = vertcat(TR{:});

lt = size(PH{1},2);
lf = size(PH{1},4);
ls = length(PH);

for r = 1:ls
    PH{r} = real(PH{r});
end


%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
PHall = {};
MCall = [];
TRall = [];

for i = 1:ls
    for r = 1:size(PH{i},1)
        PHall{count,1} = squeeze(PH{i}(r,:,:,:));       % [time x trials x freq]
        MCall(count,:) = MC{i}(r,:);
        TRall{count,1} = TR{i};
        count = count + 1;
    end
end
PH = PHall;
MC = MCall;
TR = TRall;
clear PHall MCall TRall

ls = length(PH);

%% SPLIT INTO CELL GROUPS AND TRIAL TYPES
MPH = cell(4,3);
VPH = cell(4,3);
Fp = cell(4,1);

for tt1 = 1:4                                                               % For each cell group
    k = (MC(:,1) == tt1);                                   
    if tt1 == 3
        k = (MC(:,1) == 0);         
    elseif tt1 == 4
        k = isnan(MC(:,1));
    end
    ph = PH(k);                                                             % Keep phases from that groupd only
    tr = TR(k);                                                             % Keep corresponding trials
    
    for i = 1:length(ph)                                                    % For each cell
        for tt2 = 1:2                                                       % For each first odor
            phtt = ph{i}(:,tr{i} == tt2,:);                                 % keeep those trials
            mphtt = circ_mean(phtt,[],2);                                   % Get circular mean over trials [time x 1 x lf]
            mphtt = permute(mphtt,[2,1,3]);                                 % Turn to [1 x time x lf]
            MPH{tt1,tt2} = cat(1,MPH{tt1,tt2}, mphtt);                      % Append
                   
            vphtt = circ_var(phtt,[],[],2);                                 % Repeat for circular phase variance
            vphtt = permute(vphtt,[2,1,3]);                
            VPH{tt1,tt2} = cat(1,VPH{tt1,tt2}, vphtt);
        end
    end
    
    for i = 1:length(ph)                                                    % Repeat for mean/variance across all trials
        mphtt = circ_mean(ph{i},[],2);                       
        mphtt = permute(mphtt,[2,1,3]);                
        MPH{tt1,3} = cat(1,MPH{tt1,3}, mphtt);
        
        vphtt = circ_var(ph{i},[],[],2);                  
        vphtt = permute(vphtt,[2,1,3]);                 
        VPH{tt1,3} = cat(1,VPH{tt1,3}, vphtt);
    end
end

%% GET MEAN PHASES FOR ALL CELLS AND ALL TRIALS
MPHall = nan(ls,lt,lf);
VPHall = nan(ls,lt,lf);

for i = 1:length(PH)                                                        % REPEAT FOR ALL CELLS AND ALL TRIALS
    mphtt = circ_mean(PH{i},[],2);                      
    mphtt = permute(mphtt,[2,1,3]);                 
    MPHall(i,:,:) = mphtt;
    
    vphtt = circ_var(PH{i},[],[],2);                         
    vphtt = permute(vphtt,[2,1,3]);                 
    VPHall(i,:,:) = vphtt;
end

%% SORT
for tt1 = 1:4                                                               % For each first odor
    k = (MC(:,1) == tt1);
    if tt1 == 3
        k = (MC(:,1) == 0);
    elseif tt1 == 4
        k = isnan(MC(:,1));
    end
    
    [Fp{tt1},order] = sort(MC(k,2));                                        % Sort fields
    for tt2 = 1:3
        MPH{tt1,tt2} = MPH{tt1,tt2}(order,:,:);
        VPH{tt1,tt2} = VPH{tt1,tt2}(order,:,:);
    end
end

[Fall,order] = sort(MC(:,2));                                               % Sort all fields from all cells
CTall = MC(order,1);                                                        % And keep their cell-group
MPHall = MPHall(order,:,:);
VPHall = VPHall(order,:,:);

%% COMPUTE PHASE DIFFERENCE BETWEEN THE TWO ODORS
% MFd = abs(circ_dist(MPH0{1},MPH0{2}));

%% PLOT AVERAGE PHASES OF ODOR-SPECIFIC MCELLS
% for f = 1:lf
%     figure('Name',['Odor-specific M-cells - ',freqtit{f}]);
%     for tt1 = 1:2                                                  % For each trial type
%         for tt2 = 1:2                                              % For each trial type (again)
%             subplot(2,2,(tt1-1)*2 + tt2); hold on;
%             plot_seq_rates(MPH{tt1,tt2}(:,:,f),[],time,timepoints,1); % Plot pooled sequences
%             plot(Fp{tt1},1:length(Fp{tt1}),'.k')
%         end
%     end
% end

%% PLOT MEAN PHASES AND VARIANCE OF NON-ODOR-SPECIFIC MCELLS ON ODOR1 VS ODOR2 TRIALS AND ON ALL TRIALS
for f = 1:lf
    figure('Name',['Non-odor-specific M-cells - ',freqtit{f}]);
    for tt = 1:3
        subplot(4,3,[tt 3+tt]);
        plot_seq_rates(MPH{3,tt}(:,:,f),[],time,timepoints,1);           
        plot(Fp{3},1:length(Fp{3}),'.k')
        xlim([time(1), timepoints(5)]);
        
        subplot(4,3,6+tt); hold on;
        A = MPH{3,tt}(:,:,f);
        plot_mrate_traces(circ_mean(A,[],1),circ_std(A,[],1)/sqrt(size(A,1)),time,timepoints,cols(tt,:));
        ylabel('Mean phase')

        subplot(4,3,9+tt); hold on;
        A = VPH{3,tt}(:,:,f);
        plot_mrate_traces(mean(A,1),std(A,[],1)/sqrt(size(A,1)),time,timepoints,cols(tt,:));
        ylabel('Phase variance')
    end
end

%% REPEAT FOR ALL CELLS STACKED
figure('Name','ALL cells and frequencies');
for f = 1:lf
    subplot(4,lf,[f lf+f]);
    plot_seq_rates(MPHall(:,:,f),[],time,timepoints,1);       % Plot pooled sequences
    plot(Fall(isnan(CTall)),find(isnan(CTall)),'ob')
    plot(Fall(~isnan(CTall)),find(~isnan(CTall)),'oy','Markerfacecolor','y')
    xlim([time(1), timepoints(5)]);
    
    subplot(4,lf,2*lf+f); hold on;
    plot_mrate_traces(circ_mean(MPHall(:,:,f),[],1),circ_std(MPHall(:,:,f),[],1)/sqrt(ls),time,timepoints,cols(tt,:));
    ylabel('Mean phase')
    
    subplot(4,lf,3*lf+f); hold on;
    plot_mrate_traces(mean(VPHall(:,:,f),1),std(VPHall(:,:,f),[],1)/sqrt(ls),time,timepoints,cols(tt,:));
    ylabel('Phase variance')
end

%% PLOT PHASE VARIANCE OF NON-ODOR-SPECIFIC MCELLS ON ODOR1 VS ODOR2 TRIALS AND ON ALL TRIALS
% for f = 1:lf
%     subplot(4,3,10:11); hold on;
%     plot_mrate_traces(circ_mean(MFd(:,:,f),[],1),circ_std(MFd(:,:,f),[],1)/sqrt(size(MFd,1)),time,timepoints,cols(tt,:));
%     ylabel('Phase difference between odors')
% end




