function Plot_over_trials_rate(ID,day,session,mcell_flag,positive_flag)

cols = cell_colors;
 
if nargin == 4
    positive_flag = 1;
end

%% LOAD AND SET RATES AND EXAMPLE CELLS
[asapfile,~] = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);

if positive_flag == 1
    modfile = fullfile(path,[videoname,'_Mcells.mat']);
elseif positive_flag == -1    
    modfile = fullfile(path,[videoname,'_Neg_Mcells.mat']);
end

load(modfile,'Mcells','trials','bins','timepoints','onbins');
load(Datafile,'V','B');
disp(asapfile);

onbins = find(onbins);

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
R = V.R(keep,:,:);                                                          % Keep their rates
mcells = Mcells(keep,:);                                                    % and their spiking characteristics
    
[Nr,lt,Ntr] = size(R);

%% PLOT
for c = 1:Nr                                                                % For each selected cell
    pref_od = mcells(c,1);                                                  % Keep its preferred odor
    mbin = mcells(c,2);                                                     % its field bin
    mpeak = mcells(c,3);                                                    % its peak rate
    si = mcells(c,4);                                                       % and SI
    
    Rc = squeeze(R(c,:,:))';                                                % Keep its rate over ALL trials [trials x bins]
     
    Rc = (Rc - mean(Rc(:,onbins),2)) ./ std(Rc(:,onbins),[],2);             % ZSCORE THE RATE OVER MODULATION BINS ONLY!!!
    Rc(isnan(Rc)) = 0;
    Rc(isinf(Rc)) = 0;
    
    %% COMPUTE MEAN RATES BY FIRSTODOR
    Mt = zeros(2,lt);
    St = zeros(2,lt);
    for tt = 1:2                                                             % For each first odor
        Mt(tt,:) = mean(Rc(trials == tt,:),1);                              % Keep the cell's mean firing rate
        St(tt,:) = SEM(Rc(trials == tt,:));                                 % Keep the cell's SEM
    end
    
    %% REORDER
    Rt2 = Rc(trorder,:);
    Rt3 = Rc(tresorder,:);
     
    %% GET SIGNIFICANCE
    p = ones(1,length(onbins));
    for b = 1:length(onbins)
        R1 = squeeze(Rc(trials == 1,onbins(b)));                            % Keep zscore rates over all trials 
        R2 = squeeze(Rc(trials == 2,onbins(b)));                            % across the two odors
        p(b) = ranksum(R1,R2);                                              % Compare their medians
    end
    [~,p] = fdr(p);
    
    %% PLOT ALL TRIALS
    figure('Name',[ID,' ',day,' ',videoname]);
    
    subplot(4,3,[1 4 7]);
    plot_seq_rates(Rc,[],bins,timepoints,2);           % Plot all trials (color plot)    
    for tr = 1:Ntr
        plot(licks{tr},tr*ones(size(licks{tr})),'.w');
    end
    line([1 1]*mbin,[1 Ntr],'color','w','linestyle','--','linewidth',2); % Add dashed line on its field bin
    title('All trials');
    xlim([bins(1) timepoints(5)]);
    
    subplot(4,3,10); hold on;
    plot_mrate_traces(mean(Rc,1),SEM(Rc),bins,timepoints,'k');
    
    %% PLOT EACH TRIAL TYPE
    subplot(4,3,[2 5 8]);
    plot_seq_rates(Rt2,[],bins,timepoints,2); % Plot trials split by type (color)
    for tt = 1:3
        k = sum(combos <= tt);
        line([bins(1) timepoints(7)],[1 1]*k+0.5,'color','w','linewidth',2); % Add lines between trial types
    end
    for tr = 1:Ntr
        plot(licks2{tr},tr*ones(size(licks2{tr})),'.w');
    end  
    mark_field(pref_od,mbin,trials);
    xlim([bins(1) timepoints(5)]);
    title(['SI = ',num2str(si)]);
    
    subplot(4,3,11); hold on;
    plot_mrate_traces(Mt,St,bins,timepoints,cols);
    plot(bins(onbins(p < 0.05)), max(Mt(:,onbins(p < 0.05)),[],1)+0.1,'.k');
    
    %% PLOT EACH TRIAL TYPE WITH CORRECTS ON TOP OF ERRORS
%     subplot(4,3,[3 6 9]);
%     plot_seq_rates(Rt3,[],[],bins,timepoints,2); %Plot trials split by type and outcome (color)
%     for tt = 1:3
%         k = sum(combos <= tt);
%         line([bins(1) timepoints(7)],[1 1]*k+0.5,'color','w','linewidth',2); % Add lines between trial types
%     end
%     for tt = 0.75:1:3.75
%         k = sum(combres <= tt);
%         line([bins(1) timepoints(7)],[1 1]*k+0.5,'color','c','linewidth',1.5,'linestyle','--'); % Add dashed lines between outcomes
%     end
%     for tr = 1:Ntr
%         plot(licks3{tr},tr*ones(size(licks3{tr})),'.w');
%     end
%     mark_field(pr_od,mbin,trials);    
%     xlim([bins(1) timepoints(6)]);    
%     title('Error vs Correct');
end
drawnow;

