function Plot_over_trials_DFF(ID,day,session,mcell_flag)

cols = cell_colors;

%% LOAD AND SET RATES AND EXAMPLE CELLS
[asapfile,~] = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
modfile = fullfile(path,[videoname,'_Mcells.mat']);
load(modfile,'Mcells','trials','timepoints');
load(Datafile,'B','V');
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
A = V.DespDFF(keep,:,:);                                                          % Keep their rates
mcells = Mcells(keep,:);                                                    % and their spiking characteristics

[Nr,lt,Ntr] = size(A);
    
%% PLOT
for c = 1:Nr                                                                % For each selected cell
    pref_od = mcells(c,1);                                                  % Keep its preferred odor
    mbin = mcells(c,2);                                                     % its field bin
    si = mcells(c,4);                                                       % and SI
    
    %% SHAPE 
    Dc = squeeze(A(c,:,:));                                                 % Keep its rate over ALL trials
    Dc = Dc';                                                               % turn to [trials x bins]
    
    %% SMOOTH
    Dc = sgolayfilt(Dc, 1, 101, [], 2);   
    
    %% COMPUTE MEAN BY FIRSTODOR
    MD = zeros(2,lt);
    SD = zeros(2,lt);
    for tt = 1:2                                                            % For each first odor
        MD(tt,:) = mean(Dc(trials == tt,:),1);                              % Keep the cell's mean firing rate
        SD(tt,:) = SEM(Dc(trials == tt,:));                                 % Keep the cell's SEM
    end
        
    %% GET SIGNIFICANCE
    pD = ones(1,length(ontime));
%     for f = 1:lf
%         for t = 1:length(ontime)
%             R1 = squeeze(Ac(trials == 1,ontime(t),f));                            % Keep zscore rates over all trials
%             R2 = squeeze(Ac(trials == 2,ontime(t),f));                            % across the two odors
%             pA(f,t) = ranksum(R1,R2);                                              % Compare their medians
%         end
%         [~,pA(f,:)] = fdr(pA(f,:));
%     end
    
    %% REORDER
    Dc2 = Dc(trorder,:);

    %% PLOT ALL TRIALS
    figure('Name',[ID,' ',day,' ',videoname]);
    
    subplot(4,3,[1 4 7]);
    plot_seq_rates(Dc,[],time,timepoints,0.3);           % Plot all trials (color plot)    
    for tr = 1:Ntr
        plot(licks{tr},tr*ones(size(licks{tr})),'.w');
    end
    line([1 1]*mbin,[1 Ntr],'color','w','linestyle','--','linewidth',2); % Add dashed line on its field bin
    title('All trials');
    xlim([time(1) timepoints(5)]);
    
    subplot(4,3,10); hold on;
    plot_mrate_traces(mean(Dc,1),SEM(Dc),time,timepoints,'k');
    
    %% PLOT EACH TRIAL TYPE
    subplot(4,3,[2 5 8]);
    plot_seq_rates(Dc2,[],time,timepoints,0.3); % Plot trials split by type (color)
    for tt = 1:3
        k = sum(combos <= tt);
        line([time(1) timepoints(7)],[1 1]*k+0.5,'color','w','linewidth',2); % Add lines between trial types
    end
    for tr = 1:Ntr
        plot(licks2{tr},tr*ones(size(licks2{tr})),'.w');
    end  
    mark_field(pref_od,mbin,trials);
    xlim([time(1) timepoints(5)]);
    title(['SI = ',num2str(si)]);
    
    subplot(4,3,11); hold on;
    plot_mrate_traces(MD,SD,time,timepoints,cols);
    plot(time(ontime(pD < 0.05)), max(MD(:,ontime(pD < 0.05)),[],1)+0.1,'.k');
end
drawnow;

