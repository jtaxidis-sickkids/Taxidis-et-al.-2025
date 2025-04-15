function Plot_over_trials_ampl(ID,day,session,mcell_flag)

cols = cell_colors;
freqtit = {'Delta (0.5-3Hz)';'Theta (4-10Hz)';'Beta (15-30Hz)';'Gamma (40-90Hz)'};

%% LOAD AND SET RATES AND EXAMPLE CELLS
[asapfile,~] = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
modfile = fullfile(path,[videoname,'_Mcells.mat']);
load(modfile,'Mcells','trials','timepoints');
load(Datafile,'B','F');
disp(asapfile);

Fs = 1000;
time = 0 : 1/Fs : (11*Fs-1)/Fs;
ontime = (time >= timepoints(1) & time < timepoints(3));
ontime = find(ontime);

%% ORGANIZE TRIAL SEQUENCE
combos = B.progress(:,1);                                                   % Keep the trial odor combos
[~,~,combos] = unique(combos);                                          % Order them based on combo
[combos,trorder] = sort(combos);                                        % Sort them by AA AB BA BB, and keep new trial order

%% SPLIT INTO CORRECT AND ERROR TRIALS
results = B.progress(:,2);                                                  % Keep the trial outcomes
results(results == 1 | results == 4) = +0.25;                           % Correct trials
results(results == 2 | results == 3) = -0.25;                           % Error trials
results = results(trorder);                                             % Sort outcomes according to trial order

combres = combos + results;                                             % Add outcome info to trial info (so AA/error=0.75, AA/corr=1.25,AB/error=1.75 etc)
[combres,tresorder] = sort(combres);                                    % Rearrange trials to split into error and corrects in each trial type
tresorder = trorder(tresorder);                                         % The final order is the re-ordering of the previous order

%% GET LICK TIMES
licks = get_licks(B.licks,timepoints);
licks2 = licks(trorder);
licks3 = licks(tresorder);

%% KEEP ONLY A SET OF CELLS
switch mcell_flag
    case 1
        keep = ~isnan(Mcells(:,1));     % Keep only Mcells (with and without odor preference)
    case -1
        keep = isnan(Mcells(:,1));      % Keep only non-Mcells
    case 0
        keep = ones(size(Mcells,1),1);  % Keep all cells
end
A = F.Ampl(keep,:,:,:);                                                          % Keep their rates
mcells = Mcells(keep,:);                                                    % and their spiking characteristics

[Nr,lt,Ntr,lf] = size(A);
    
%% PLOT
for c = 1:Nr                                                         % For each selected cell
    pref_od = mcells(c,1);                                                  % Keep its preferred odor
    mbin = mcells(c,2);                                                     % its field bin
    si = mcells(c,4);                                                       % and SI
    
    %% SHAPE AMPLITUDE AND PHASES
    Ac = squeeze(A(c,:,:,:));                                               % Keep its rate over ALL trials
    Ac = permute(Ac,[2 1 3]);                                               % turn to [trials x bins x freqs]
    
    %% ZSCORE AMPLITUDE
    Ac = (Ac - mean(Ac(:,ontime,:),2)) ./ std(Ac(:,ontime,:),[],2);         % ZSCORE OVER MODULATION BINS ONLY!!!
    Ac(isnan(Ac)) = 0;
    Ac(isinf(Ac)) = 0;
    
    %% COMPUTE MEAN BY FIRSTODOR
    MA = zeros(2,lt,lf);
    SA = zeros(2,lt,lf);
    for tt = 1:2                                                             % For each first odor
        MA(tt,:,:) = mean(Ac(trials == tt,:,:),1);                              % Keep the cell's mean firing rate
        for f = 1:lf
            SA(tt,:,f) = SEM(Ac(trials == tt,:,f));                                 % Keep the cell's SEM
        end
    end
        
    %% GET SIGNIFICANCE
    pA = ones(lf,length(ontime));
%     for f = 1:lf
%         for t = 1:length(ontime)
%             R1 = squeeze(Ac(trials == 1,ontime(t),f));                            % Keep zscore rates over all trials
%             R2 = squeeze(Ac(trials == 2,ontime(t),f));                            % across the two odors
%             pA(f,t) = ranksum(R1,R2);                                              % Compare their medians
%         end
%         [~,pA(f,:)] = fdr(pA(f,:));
%     end
    
    %% REORDER
    Ac2 = Ac(trorder,:,:);
    
    %% PLOT AMPLITUDES
    figure('Name',[ID,' ',day,' ',videoname]); 
    for f = 1:lf
        %% PLOT OVER ALL TRIALS
        subplot(4,8,(f-1)*2+[1 9 17]);
        plot_seq_rates(Ac(:,:,f),[],time,timepoints,2);           % Plot all trials (color plot)
        for tr = 1:Ntr
            plot(licks{tr},tr*ones(size(licks{tr})),'.w');
        end
        line([1 1]*mbin,[1 Ntr],'color','w','linestyle','--','linewidth',2); % Add dashed line on its field bin
        title(freqtit{f});
        xlim([0 timepoints(5)]);
        
        subplot(4,8,(f-1)*2+25); hold on;
        plot_mrate_traces(mean(Ac(:,:,f),1),SEM(Ac(:,:,f)),time,timepoints,'k');
        
        %% PLOT EACH TRIAL TYPE
        subplot(4,8,(f-1)*2+[2 10 18]);
        plot_seq_rates(Ac2(:,:,f),[],time,timepoints,2); % Plot trials split by type (color)
        for tt = 1:3
            k = sum(combos <= tt);
            line([0 timepoints(7)],[1 1]*k+0.5,'color','w','linewidth',2); % Add lines between trial types
        end
        for tr = 1:Ntr
            plot(licks2{tr},tr*ones(size(licks2{tr})),'.w');
        end
        mark_field(pref_od,mbin,trials);
        xlim([0 timepoints(5)]);
        
        subplot(4,8,(f-1)*2+26); hold on;
        plot_mrate_traces(MA(:,:,f),SA(:,:,f),time,timepoints,cols);
        plot(time(ontime(pA < 0.05)), max(MA(:,ontime(pA < 0.05),f),[],1)+0.1,'.k');
    end
end
drawnow;

