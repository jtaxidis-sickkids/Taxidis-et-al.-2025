function Theta_Rates_Locomotion(INclass)

%% LOAD ALL MOTION SEGMENTS IN ALL FILES
if strcmp(INclass,'PV')
    Days_PV;
    maxd = 7;
elseif strcmp(INclass,'SOM')
    Days_SOM;
    maxd = 8;
end

Na = length(asap);
MOV = cell(Na,maxd);
MOV_t = cell(Na,maxd);

for a = 1:Na
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1:length(days)                                               % For each selected day      
        MOV{a,d} = {};
        MOV_t{a,d} = {};

        for s = 1:length(ses{d})                                    % for each selected session that day
            [asapfile,~] = get_ASAPfile(ID,days{d},ses{d}(s));
            [path,videoname] = fileparts(asapfile);
            disp(asapfile);
 
            datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
            load(datafile,'M');           
            
            MOV{a,d} = cat(1,MOV{a,d},M.Move);            
            MOV_t{a,d} = cat(1,MOV_t{a,d},{M.Move_time'});            
        end
    end
end

%% LOAD POOLED DATA
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'AMP','R','DFF','MOT','bins','time','timepoints');

%% GET EDGES AND ODOR TIMES AND BINS
t_no_edges = (time > 0.7 & time < 10.3);
b_no_edges = (bins > 0.7 & bins < 10.3);

od1 = (time >= timepoints(1) & time <= timepoints(2));
od2 = (time >= timepoints(3) & time <= timepoints(4));
t_no_od = ~(od1 | od2);

od1 = (bins >= timepoints(1) & bins <= timepoints(2));
od2 = (bins >= timepoints(3) & bins <= timepoints(4));
b_no_od = ~(od1 | od2);

%% KEEP ONLY NAIVE OR TRAINED DAYS
% AMP = AMP(:,1:2);
% R = R(:,1:2);
% MOT = MOT(:,1:2);
% MOV = MOV(:,1:2);
% MOV_t = MOV_t(:,1:2);

% AMP(:,1:2) = [];
% R(:,1:2) = [];
% MOT(:,1:2) = [];
% MOV(:,1:2) = [];
% MOV_t(:,1:2) = [];

%% STACK IN SINGLE COLUMNS
R = vertcat(R{:});
DFF = vertcat(DFF{:});
AMP = vertcat(AMP{:});
MOT = vertcat(MOT{:});
MOV = vertcat(MOV{:});
MOV_t = vertcat(MOV_t{:});

ls = length(R); % Number of sessions

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
Rall = [];
Dall = [];
AMPall = [];
MOTall = {};
MOVall = {};
MOV_t_all = {};

for i = 1:ls
    for r = 1:size(R{i},1)
        Rall{count,1} = squeeze(R{i}(r,:,:))';                               % Keep only that cell's rates [trials x bins]
        Dall{count,1} = squeeze(DFF{i}(r,:,:))';                               % Keep only that cell's DFF [trials x bins]
        AMPall{count,1} = squeeze(AMP{i}(r,:,:,2))';                         % Keep only that cell's rates THETA amplitude [trials x time]
        MOTall{count,1} = MOT{i}';                                           % (Same dimensions)
        MOVall{count,1} = MOV{i};                                           % (Same dimensions)
        MOV_t_all{count,1} = MOV_t{i};                                           % (Same dimensions)
        count = count + 1;
    end
end
R = Rall;
DFF = Dall;
AMP = AMPall;
MOT = MOTall;
MOV = MOVall;
MOV_t = MOV_t_all;
clear Rall Dall MOTall MOVall MOV_t_all AMPall

ls = count-1;
Ntr = cellfun(@length, MOV_t); % Number of sessions

%% SMOOTH DFF RATES AND MOTION
for i = 1:ls
    R{i} = smoothdata(R{i},2,'movmean',5);
    MOT{i} = smoothdata(MOT{i},2,'movmean',50);
    DFF{i} = sgolayfilt(DFF{i}, 1, 21, [] , 2);                             % Apply linear filter along time bins
end

%% ZSCORE THETA, RATES AND LOCOMOTION
zAMP = cell(ls,1);
zMOT = cell(ls,1);
zR = cell(ls,1);
zD = cell(ls,1);

for i = 1:ls
    zAMP{i} = zscore(AMP{i},[],2);
    zMOT{i} = zscore(MOT{i},[],2);
    zR{i} = zscore(R{i},[],2);
    zD{i} = zscore(DFF{i},[],2);
end

%% GET MEANS ACROSS TRIALS 
mzAMP = cellfun(@(x) mean(x,1),zAMP,'UniformOutput',false);
mzAMP = cell2mat(mzAMP);

mzMOT = cellfun(@(x) mean(x,1),zMOT,'UniformOutput',false);
mzMOT = cell2mat(mzMOT);

mzR = cellfun(@(x) mean(x,1),zR,'UniformOutput',false);
mzR = cell2mat(mzR);

mzD = cellfun(@(x) mean(x,1),zD,'UniformOutput',false);
mzD = cell2mat(mzD);

%% PLOT MEANS ACROSS TRIAL
figure;
fill_plot(time,mean(mzAMP,1),SEM(mzAMP),'k')
hold on;
fill_plot(time,mean(mzMOT,1),SEM(mzMOT),'r')
fill_plot(time,mean(mzD,1),SEM(mzD),'m')
fill_plot(bins,mean(mzR,1),SEM(mzR),'b')
xlim([0.7 10.3]);

%% FIND MOTION BINS
lb = length(bins);

MOV_b = cell(ls,1);
for i = 1:ls
    MOV_b{i} = zeros(Ntr(i),lb);
end

for i = 1:ls
    for tr = 1:Ntr(i)
        mot = [];
        for m = 1:size(MOV_t{i}{tr},1)                                                 % For each motion segment
            [~,k1] = min(abs(bins - MOV_t{i}{tr}(m,1)));                               % Find the bin closest to motion initiation
            [~,k2] = min(abs(bins - MOV_t{i}{tr}(m,2)));                               % and termination
            k2 = min(k2,lb);                                % In case it goes just above lb
            mot = [mot, k1:k2];                                                     % Keep the motion indexes
        end
        MOV_b{i}(tr,mot) = 1;
    end
end

%% SPLIT INTO LOCOMOTION AND IMMOBILITY SEGMENTS
mAMPm = cell(ls,1);
mAMPi = cell(ls,1);
mRm = cell(ls,1);
mRi = cell(ls,1);

for i = 1:ls
    for tr = 1:Ntr(i)
        k = MOV{i}(tr,:);
        k1 = k & t_no_edges;
        k2 = ~k & t_no_edges;
        
        mAMPm{i}(tr) = mean(AMP{i}(tr,k1));
        mAMPi{i}(tr) = mean(AMP{i}(tr,k2));
        
        k = MOV_b{i}(tr,:);
        k1 = k & b_no_edges;
        k2 = ~k & b_no_edges;
        
        mRm{i}(tr) = mean(R{i}(tr,k1));
        mRi{i}(tr) = mean(R{i}(tr,k2));
    end
end
mAMPm = cell2mat(mAMPm');
mAMPi = cell2mat(mAMPi'); 
mRm = cell2mat(mRm');
mRi = cell2mat(mRi');

%% COMPUTE MOTIONvsIMMOBILITY DFF SPECTRA
Sm = nan(ls,8001);
Si = nan(ls,8001);                                                          % 4001 = length of Freq vector

for i = 1:ls
    Dconcat = DFF{i}';            % turn DFF to [time x trials]
    Dconcat = Dconcat(:);            % Concatenate all trials in one column

    no_ed = repmat(t_no_edges,[Ntr(i),1]);
    no_ed = no_ed';
    no_ed = no_ed(:);
    
    k = MOV{i}';
    k = k(:);
    
    k1 = k & no_ed;
    k2 = ~k & no_ed;

    if sum(k1) >= 5000      % IF THERE ARE AT LEAST 5 sec OF MOTION
        [Sm(i,:),Fvec] = pspectrum(Dconcat(k1)',1000,'FrequencyLimits',[0 500],'FrequencyResolution',0.5);      % Compute power spectrum during motion
    end
    if sum(k2) >= 5000      % IF THERE ARE AT LEAST 5 sec OF IMMOBILITY
        [Si(i,:),Fvec] = pspectrum(Dconcat(k2)',1000,'FrequencyLimits',[0 500],'FrequencyResolution',0.5);      % Compute power spectrum during motion
    end
end

Sm(~isnan(Sm)) = pow2db(Sm(~isnan(Sm)));                                    % Switch to decibel (skip Nan entries where it crashes)
Si(~isnan(Si)) = pow2db(Si(~isnan(Si)));

%% COMPUTE NONODOR MOTIONvsIMMOBILITY DFF SPECTRA
Smno = nan(ls,8001);
Sino = nan(ls,8001);                                                          % 4001 = length of Freq vector

for i = 1:ls
    Dconcat = DFF{i}';            % turn DFF to [time x trials]
    Dconcat = Dconcat(:);            % Concatenate all trials in one column

    no_ed = repmat(t_no_edges,[Ntr(i),1]);
    no_ed = no_ed';
    no_ed = no_ed(:);
    
    no_od = repmat(t_no_od,[Ntr(i),1]);
    no_od = no_od';
    no_od = no_od(:);
    
    k = MOV{i}';
    k = k(:);
    
    k1 = k & no_ed & no_od;
    k2 = ~k & no_ed & no_od;

    if sum(k1) >= 5000      % IF THERE ARE AT LEAST 5 sec OF MOTION
        [Smno(i,:),Fvec] = pspectrum(Dconcat(k1)',1000,'FrequencyLimits',[0 500],'FrequencyResolution',0.5);      % Compute power spectrum during motion
    end
    if sum(k2) >= 5000      % IF THERE ARE AT LEAST 5 sec OF IMMOBILITY
        [Sino(i,:),Fvec] = pspectrum(Dconcat(k2)',1000,'FrequencyLimits',[0 500],'FrequencyResolution',0.5);      % Compute power spectrum during motion
    end
end

Smno(~isnan(Smno)) = pow2db(Smno(~isnan(Smno)));                                    % Switch to decibel (skip Nan entries where it crashes)
Sino(~isnan(Sino)) = pow2db(Sino(~isnan(Sino)));

%% PLOT AVERAGE SPECTRA
figure;
subplot(121); hold on;
fill_plot(Fvec',nanmean(Sm,1),SEM(Sm),'r'); 
fill_plot(Fvec',nanmean(Si,1),SEM(Si),'k'); 
ylabel('Spectrum motion vs immobility');
set(gca, 'YScale', 'log')

% Significance per frequency
pS = ones(size(Fvec));
for i = 1:length(Fvec)
    pS(i) = ranksum(Sm(:,i),Si(:,i));
end
[~,pS] = fdr(pS);
plot(Fvec(pS < 0.05),0.01*ones(1,sum(pS<0.05)),'*k')

subplot(122); hold on;
fill_plot(Fvec',nanmean(Smno,1),SEM(Smno),'r'); 
fill_plot(Fvec',nanmean(Sino,1),SEM(Sino),'k'); 
ylabel('Spectrum motion vs immobility excluding odors');
set(gca, 'YScale', 'log')

% Significance per frequency
pS = ones(size(Fvec));
for i = 1:length(Fvec)
    pS(i) = ranksum(Smno(:,i),Sino(:,i));
end
[~,pS] = fdr(pS);
plot(Fvec(pS < 0.05),0.01*ones(1,sum(pS<0.05)),'*k')


%% COMPARE MOTION vs IMMOB THETA & RATES
[~,p] = ttest(mAMPi,mAMPm);           % Paired sample ttest
[~,pR] = ttest(mRi,mRm);           % Paired sample ttest

figure;
subplot(321)
scatter_plot_comparison(mAMPi,mAMPm,p);
xlabel('theta imm');
ylabel('theta mot');
title(num2str(p));

subplot(322)
scatter_plot_comparison(mRi,mRm,pR);
xlabel('rates imm');
ylabel('rates mot');
title(num2str(pR));

%% SPLIT INTO LOCOMOTION AND IMMOBILITY BUT EXCLUDING ODOR SEGMENTS
mAMPm = cell(ls,1);
mAMPi = cell(ls,1);
mRm = cell(ls,1);
mRi = cell(ls,1);

for i = 1:ls
    for tr = 1:Ntr(i)
        k = MOV{i}(tr,:) ;
        k1 = k & t_no_od & t_no_edges;
        k2 = ~k & t_no_od & t_no_edges;
        
        mAMPm{i}(tr) = mean(AMP{i}(tr,k1));
        mAMPi{i}(tr) = mean(AMP{i}(tr,k2));
        
        k = MOV_b{i}(tr,:);
        k1 = k & b_no_od & b_no_edges;
        k2 = ~k & b_no_od & b_no_edges;
        
        mRm{i}(tr) = mean(R{i}(tr,k1));
        mRi{i}(tr) = mean(R{i}(tr,k2));
    end
end
mAMPm = cell2mat(mAMPm');
mAMPi = cell2mat(mAMPi');
mRm = cell2mat(mRm');
mRi = cell2mat(mRi');

%% COMPARE MOTION vs IMMOB THETA & RATES
[~,p] = ttest(mAMPi,mAMPm);           % Paired sample ttest
[~,pR] = ttest(mRi,mRm);           % Paired sample ttest

subplot(323)
scatter_plot_comparison(mAMPi,mAMPm,p);
xlabel('theta imm no odor');
ylabel('theta mot no odor');
title(num2str(p));

subplot(324)
scatter_plot_comparison(mRi,mRm,pR);
xlabel('rates imm no odor');
ylabel('rates mot no odor');
title(num2str(pR));

%% COMPARE ODOR vs NO-ODOR THETA & RATES
mAMPo = cell(ls,1);
mAMPno = cell(ls,1);
mMOTo = cell(ls,1);
mMOTno = cell(ls,1);
        
k1 = t_no_od & t_no_edges;
k2 = ~t_no_od;
        
for i = 1:ls
    for tr = 1:size(R{i},1)
        mAMPno{i}(tr) = mean(AMP{i}(tr,k1));
        mAMPo{i}(tr) = mean(AMP{i}(tr,k2));
            
        mMOTno{i}(tr) = mean(MOT{i}(tr,k1));
        mMOTo{i}(tr) = mean(MOT{i}(tr,k2));
    end
end

mAMPo = cell2mat(mAMPo');
mAMPno = cell2mat(mAMPno');
mMOTo = cell2mat(mMOTo');
mMOTno = cell2mat(mMOTno');

%% COMPARE MOTION vs IMMOB THETA 
[~,p] = ttest(mAMPno,mAMPo);           % Paired sample ttest
[~,pM] = ttest(mMOTno,mMOTo);           % Paired sample ttest

subplot(325)
scatter_plot_comparison(mAMPno,mAMPo,p);
xlabel('theta no odor');
ylabel('theta odor');
title(num2str(p));

subplot(326)
scatter_plot_comparison(mMOTno,mMOTo,pM);
xlabel('motion no odor');
ylabel('motion odor');
title(num2str(pM));
