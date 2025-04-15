function Plot_over_trials_phase(ID,day,session,mcell_flag)

cols = cell_colors;
freqtit = {'Delta (1-4Hz)';'Theta (5-11Hz)';'Beta (15-30Hz)';'Gamma (40-90Hz)'};
addpath('CircStat2012a');

%% LOAD AND SET RATES AND EXAMPLE CELLS
[asapfile,~] = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
modfile = fullfile(path,[videoname,'_Mcells.mat']);
load(modfile,'Mcells','trials','timepoints');
load(Datafile,'B','F','V');
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

P = F.Phase(keep,:,:,:);
S = V.SP(keep,:);
SP = F.SP_ph(keep,:);
mcells = Mcells(keep,:);                                                    % and their spiking characteristics

[Nr,lt,Ntr,lf] = size(P);

%% PLOT
for c = 1:Nr                                                         % For each selected cell
    pref_od = mcells(c,1);                                                  % Keep its preferred odor
    mbin = mcells(c,2);                                                     % its field bin
    
    %% SHAPE PHASES
    Pc = squeeze(P(c,:,:,:));                                               % Keep its rate over ALL trials
    Pc = permute(Pc,[2 1 3]);                                               % turn to [trials x bins x freqs]
    
    Sc = S(c,:);
    SPc = SP(c,:);
    
    %% COMPUTE MEAN BY FIRSTODOR
    MPh = zeros(3,lt,lf);
    SPh = zeros(3,lt,lf);
    for f = 1:lf
        for tt = 1:2                                                            % For each first odor
            MPh(tt,:,f) = circ_mean(Pc(trials == tt,:,f),[],1);                  % Keep the cell's mean firing rate
            SPh(tt,:,f) = circ_std(Pc(trials == tt,:,f),[],[],1) / sqrt(sum(trials == tt));            % Keep the cell's SEM
        end
        
        % ADD PHASE DISFFERENCE BETWEEN THE TWO ODORS
        MPh(3,:,f) = abs(circ_dist(MPh(1,:,f),MPh(2,:,f)));  
    end
    
    
    %% COMPUTE MEAN OVER ALL TRIALS
    MPH = circ_mean(Pc,[],1);
    SPH = circ_std(Pc,[],[],1) / sqrt(size(Pc,1));            % Keep the cell's SEM
    
    %% GET SIGNIFICANCE
    pP = ones(lf,length(ontime));
    %     for f = 1:lf
    %         for t = 1:length(ontime)
    %             R1 = squeeze(Pc(trials == 1,ontime(t),f));                            % Keep zscore rates over all trials
    %             R2 = squeeze(Pc(trials == 2,ontime(t),f));                            % across the two odors
    %             pP(f,t) = ranksum(R1,R2);                                              % Compare their medians
    %         end
    %         [~,pP(f,:)] = fdr(pP(f,:));
    %     end
    
    %% REORDER
    Pc2 = Pc(trorder,:,:);
    
    %% PLOT PHASES
    figure('Name',[ID,' ',day,' ',videoname]);
    %% PLOT OVER ALL TRIALS
    for f = 1:lf
        subplot(8,2,(f-1)*4+1);
        plot_seq_rates(Pc(:,:,f),[],time,timepoints,2);           % Plot all trials (color plot)
        for tr = 1:Ntr
            plot(licks{tr},tr*ones(size(licks{tr})),'.w');
        end
        line([1 1]*mbin,[1 Ntr],'color','w','linestyle','--','linewidth',2); % Add dashed line on its field bin
        title(freqtit{f});
        xlim([0 timepoints(5)]);
        
        subplot(8,2,(f-1)*4+2); hold on;
        plot_mrate_traces(MPH(:,:,f),SPH(:,:,f),time,timepoints,'k');
    end
    
    %% PLOT EACH TRIAL TYPE
    for f = 1:lf
        subplot(8,2,(f-1)*4+3);
        plot_seq_rates(Pc2(:,:,f),[],time,timepoints,2); % Plot trials split by type (color)
        for tt = 1:3
            k = sum(combos <= tt);
            line([0 timepoints(7)],[1 1]*k+0.5,'color','w','linewidth',2); % Add lines between trial types
        end
        for tr = 1:Ntr
            plot(licks2{tr},tr*ones(size(licks2{tr})),'.w');
        end
        mark_field(pref_od,mbin,trials);
        xlim([0 timepoints(5)]);
        
        subplot(8,2,f*4); hold on;
        plot(time,Pc(:,:,f),'Color',[.8 .8 .8]);          % SUPERIMPOSE ALL TRIALS
%         plot_mrate_traces(MPh(:,:,f),SPh(:,:,f),time,timepoints,cols);
        plot(time(ontime(pP < 0.05)), max(MPh(:,ontime(pP < 0.05),f),[],1)+0.1,'.k');
        xlim([0 timepoints(5)]);  
    end
    
    %% PLOT SPIKE-PHASE DISTRIBUTION
    figure('Name',[ID,' ',day,' ',videoname]);
    Cr = rand(Ntr,3);
    for f = 1:lf
        subplot(lf,5,(f-1)*5 + (1:2)); hold on;
        for tr = 1:Ntr
            plot(Sc{tr},SPc{tr}(:,f),'.','Color',Cr(tr,:))
        end
        axis tight;

        spc_all = cell2mat(SPc');      
        sc_all = cell2mat(Sc)';      
        
        subplot(lf,5,(f-1)*5 + (3:4)); hold on;
        t = 0:0.5:11;
        ph = -pi:0.1:pi;
        h = hist3([sc_all, spc_all(:,f)],'Edges',{t, ph});
        h = h';
        
        % Interpolate
        tq = t(1):0.1:t(end);
        phq = ph(1):0.01:ph(end);
        [T,Fint] = meshgrid(t,ph);
        [Tq,Fq] = meshgrid(tq,phq);
        h = interp2(T,Fint,h,Tq,Fq);
        imagesc(0:0.2:11, -pi:0.2:pi, h);
        axis tight;
        
        subplot(lf,5,(f-1)*5+5);hold on; view(90, -90)
        histogram(spc_all(:,f), -pi:0.2:pi);   
    end
    
end
drawnow;

