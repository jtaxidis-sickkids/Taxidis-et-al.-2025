function PERF = Figure_Progress(INclass,cell_flag)  

if nargin == 1, cell_flag = 'all'; end

%% LOAD POOLED DATA
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'R','MC','AMP','MDIR','LVEC','MOT','PROG','time','bins');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    R{1,5} = [];    R{2,3} = [];
    MC{1,5} = [];   MC{2,3} = [];
    AMP{1,5} = [];  AMP{2,3} = [];
    MDIR{1,5} = []; MDIR{2,3} = [];
    LVEC{1,5} = []; LVEC{2,3} = [];
    MOT{1,5} = [];  MOT{2,3} = [];
    PROG{1,5} = []; PROG{2,3} = [];
end

% REMOVE EXRA NAIVE DAY FROM ASAP3-1 TO GET 2 NAIVE DAYS FOR ALL MICE 
% AND SHIFT THE TRAINED-DAYS BACK SO THERE IS NO GAP-DAY
if strcmp(INclass,'SOM')
    R(3,3:7) = R(3,4:8);
    MC(3,3:7) = MC(3,4:8);
    AMP(3,3:7) = AMP(3,4:8);
    MDIR(3,3:7) = MDIR(3,4:8);
    LVEC(3,3:7) = LVEC(3,4:8);
    MOT(3,3:7) = MOT(3,4:8);
    PROG(3,3:7) = PROG(3,4:8);
end

%% INITIALIZE
[Na,maxd] = size(R);
lf = 4;

LMC = nan(Na,maxd);
NR = nan(Na,maxd);
SI = cell(Na,maxd);
PERF = cell(Na,maxd);
REJ = cell(Na,maxd);
MR = cell(Na,maxd);
MF = cell(Na,maxd);
MD = cell(Na,maxd);
LV = cell(Na,maxd);

%% LOAD AND ORGANIZE DATA
for a = 1:Na
    for d = 1:maxd
        if ~isempty(R{a,d})
            mcells = cell2mat(MC{a,d});
            if isempty(mcells), mcells = zeros(0,4); end
            
            NR(a,d) = size(mcells,1);                                       % Number of all cells
            MOT{a,d} = cellfun(@(x) mean(x,'all'), MOT{a,d});               % Get mean locomotion per session
            
            % CHOOSE CELL GROUP TO USE
            if strcmp(cell_flag,'all')
                k = true(NR(a,d),1);
            elseif strcmp(cell_flag,'mcells')
                k = ~isnan(mcells(:,1));
            elseif strcmp(cell_flag,'odorcells')
               k = mcells(:,1) > 0; 
            elseif strcmp(cell_flag,'nonodorcells')
                k = mcells(:,1) == 0;
            elseif strcmp(cell_flag,'nonmcells')
                k = isnan(mcells(:,1));
            end
            
            % KEEP CELL GROUP DATA
            LMC(a,d) = sum(k);                                              % Count all cells of the selected group
            SI{a,d} = abs(mcells(k,4));                                     % SI of selected cells
            
%             ontime = time >=1 & time <= 8; 
%             onbins = bins >=1 & bins <= 8;
%             for i = 1:size(R{a,d},1)
%                 R{a,d}{i} = R{a,d}{i}(:,onbins,:);
%                 AMP{a,d}{i} = AMP{a,d}{i}(:,ontime,:,:);
%             end
            
            % GET MEAN RATE, PHASES AND AMPLITUDES OF ALL CELLS
            MR{a,d} = cell2mat(cellfun(@(x) mean(x,[2,3]), R{a,d}, 'UniformOutput', false)); % Mean rate across all bins and trials for each cell
            MD{a,d} = cell2mat(MDIR{a,d});                                  % Combine preferred phases of all videos of that day
            LV{a,d} = cell2mat(LVEC{a,d});                                  % Same for stength of phase locking 
            MF{a,d} = cell2mat(cellfun(@(x) mean(x,[2,3]), AMP{a,d}, 'UniformOutput', false)); % Mean amplitude across all bins and trials for each cell and frequency
            MF{a,d} = reshape(MF{a,d},[NR(a,d),lf]);                        % Set to [cells x frequencies])
            
            % KEEP ONLY SELECTED CELL GROUP
            MR{a,d} = MR{a,d}(k);  
            MD{a,d} = MD{a,d}(k,:);                                         
            LV{a,d} = LV{a,d}(k,:);
            MF{a,d} = MF{a,d}(k,:);                                        
            
            % GET PERFORMANCE AND REJECTION RATES
            for s = 1:length(PROG{a,d})
                prog = PROG{a,d}{s};
                Ntr = size(prog,1);
                
                combos = prog(:,1);
                match = (rem(combos,11) == 0);
                outcomes = prog(:,2);
                PERF{a,d}(s,1) = sum(outcomes == 1 | outcomes == 4) / Ntr * 100;
                REJ{a,d}(s,1) = sum(outcomes == 4) / sum(match) * 100;
            end
        end
    end
end

LMC = LMC ./ NR *100;                       % Turn to % Mcells

% SET ALL PERFORMANCES AT NAIVE STAGE TO 50 (AND REJECTIONS TO 100
% TO AVOID WEIRD NUMBERS DUE TO TRIAL NUMBER (AND THE STEP2 THAT WAS GIVEN TO ONE SOM MOUSE)
for a = 1:Na
    for d = 1:2
        PERF{a,d} = 50 * ones(size(PERF{a,d}));
        REJ{a,d} = 100 * ones(size(REJ{a,d}));
    end
end

S = {PERF; MOT; MR; SI};
Slabel = {'Performance';'Motion';'Rates';'SI'};
freqtit = {'Delta (1-4Hz)';'Theta (5-11Hz)';'Beta (15-30Hz)';'Gamma (40-90Hz)'};
ls = length(S);
cols = 12;

%% PLOT PERFORMANCE,LOCOMOTION,RATES,SI ACROSS DAYS AND COMPARE NAIVE-TRAINED
figure('Name',[INclass,'  ',cell_flag]);

for i = 1:ls                                                                % For each measure
    % CONCATENATE MICE
    A = cell(1,maxd);
    for d = 1:maxd                                                          % For each day
        atemp = S{i}(:,d);                                                  % Keep all mice
        atemp = cell2mat(atemp);                                            % Turn to matrix 
        A{d} = atemp;                                                       % Store as cell entry
    end

    k = cellfun(@length,A);                                                 % Count #cells in each day
    kmax = max(k);                                                          % Max # cells on a single day
    D = nan(kmax,maxd);
    for d = 1:maxd                                  
        D(1:k(d),d) = A{d};
    end
    
    [Rs,Ps] = Spearman_over_days(D);
    
    subplot(ls+1,cols,cols*(i-1)+(1:2)); hold on;
    for d = 1:maxd
        plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , D(:,d),'o','Color',[0.8 0.8 0.8]);
    end
    fill_plot(1:maxd,nanmean(D),SEM(D),'b')
    title([num2str(Rs),' (',num2str(Ps),')']);
    ylabel(Slabel{i});
    xlim([0.5 maxd+0.5]);
    
    subplot(ls+1,cols,cols*(i-1)+3);
    plot_naive_trained(D);
end

%% PLOT NUMBER OF MCELLS ACROSS DAYS AND COMPARE NAIVE-TRAINED
[Rs,Ps] = Spearman_over_days(LMC);

subplot(ls+1,cols,cols*ls+(1:2));
plot_over_days(LMC);
ylabel('%Mcells');
title([num2str(Rs),' (',num2str(Ps),')']);
xlim([0.5 maxd+0.5]);
    
subplot(ls+1,cols,cols*ls+3);
plot_naive_trained(LMC);

%% PLOT BANDPOWER ACROSS DAYS AND COMPARE NAIVE-TRAINED
for f = 1:lf                                                                % For each frequency
    % CONCATENATE MICE
    A = cell(1,maxd);
    for d = 1:maxd                                                          % For each day
        atemp = MF(:,d);                                                    % Keep all mice
        atemp = cell2mat(atemp);                                            % Turn to matrix 
        A{d} = atemp(:,f);                                                  % Store are cell entry
    end

    k = cellfun(@length,A);                                                 % Count #cells in each day
    kmax = max(k);                                                          % Max # cells on a single day
    D = nan(kmax,maxd);
    for d = 1:maxd                                  
        D(1:k(d),d) = A{d};
    end
    
    [Rs,Ps] = Spearman_over_days(D);
    
    subplot(ls+1,cols,cols*(f-1)+(4:5)); hold on;
    for d = 1:maxd
        plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , D(:,d),'o','Color',[0.8 0.8 0.8]);
    end
    fill_plot(1:maxd,nanmean(D),SEM(D),'b')
    title([num2str(Rs),' (',num2str(Ps),')']);
    ylabel(freqtit{f});
    xlim([0.5 maxd+0.5]);
        
    subplot(ls+1,cols,cols*(f-1)+6);
    plot_naive_trained(D);
end

%% PLOT STRENGTH OF PHASE LOCKING ACROSS DAYS AND COMPARE NAIVE-TRAINED
for f = 1:lf                                                                % For each frequency
    % CONCATENATE MICE
    A = cell(1,maxd);
    for d = 1:maxd                                                          % For each day
        atemp = LV(:,d);                                                    % Keep all mice
        atemp = cell2mat(atemp);                                            % Turn to matrix 
        A{d} = atemp(:,f);                                                  % Store are cell entry
    end

    k = cellfun(@length,A);                                                 % Count #cells in each day
    kmax = max(k);                                                          % Max # cells on a single day
    D = nan(kmax,maxd);
    for d = 1:maxd                                  
        D(1:k(d),d) = A{d};
    end
    
    [Rs,Ps] = Spearman_over_days(D);
    
    subplot(ls+1,cols,cols*(f-1)+(7:8)); hold on;
    for d = 1:maxd
        plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , D(:,d),'o','Color',[0.8 0.8 0.8]);
    end
    fill_plot(1:maxd,nanmean(D),SEM(D),'b')
    title([num2str(Rs),' (',num2str(Ps),')']);
    ylabel(freqtit{f});
    xlim([0.5 maxd+0.5]);
        
    subplot(ls+1,cols,cols*(f-1)+9);
    plot_naive_trained(D);
end

%% PLOT PREFERRED PHASE ACROSS DAYS AND COMPARE NAIVE-TRAINED
for f = 1:lf                                                                % For each frequency
    % CONCATENATE MICE
    A = cell(1,maxd);
    for d = 1:maxd                                                          % For each day
        atemp = MD(:,d);                                                    % Keep all mice
        atemp = cell2mat(atemp);                                            % Turn to matrix 
        A{d} = atemp(:,f);                                                  % Store as cell entry
    end

    k = cellfun(@length,A);                                                 % Count #cells in each day
    kmax = max(k);                                                          % Max # cells on a single day
    D = nan(kmax,maxd);
    md = nan(1,maxd);
    sd = nan(1,maxd);
    for d = 1:maxd   
        D(1:k(d),d) = A{d};
        md(d) = circ_mean(A{d});
        sd(d) = circ_std(A{d})/sqrt(length(A{d}));
    end
    
    [Rs,Ps] = Spearman_over_days(D);
    
    subplot(ls+1,cols,cols*(f-1)+(10:11)); hold on;
    for d = 1:maxd
        plot(d*ones(kmax,1) + randn(kmax,1)*0.05 , D(:,d),'o','Color',[0.8 0.8 0.8]);
    end
    fill_plot(1:maxd,md,sd,'b');     
    title([num2str(Rs),' (',num2str(Ps),')']);
    ylabel(freqtit{f});
    xlim([0.5 maxd+0.5]);
        
    Dpre = D(:,1:2);
    Dpost = D(:,3:end);
    
    subplot(ls+1,cols,cols*(f-1)+12);
    dp = table(Dpre,Dpost,'VariableNames',{'Naive','Trained'});
    violinplot(dp);                                             % REPLACE WITH ACTUAL MEDIANS ETC
    hold on;
    
    Dpre = Dpre(:);
    Dpost = Dpost(:);
    Dpre(isnan(Dpre)) = [];
    Dpost(isnan(Dpost)) = [];
    p_prepost = circ_wwtest(Dpre,Dpost);                          % Compare their mean phases with parametric Watson-Williams test
    title(num2str(p_prepost));
end

drawnow;
