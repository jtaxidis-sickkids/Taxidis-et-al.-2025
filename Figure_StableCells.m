function Figure_StableCells(INclass)

if strcmp(INclass,'PV')
    Days_PV
    maxd = 7;
elseif strcmp(INclass,'SOM')
    Days_SOM
    maxd = 8;
end

Na = length(asap);                                                          % total number of animals
Nr = length([asap.samecells]');                                             % total number of cells followed across days

C = cell_colors;

count = 1;

%% PLOT RATES AND FOV OF REGISTERED CELLS ACROSS DAYS
for a = 1:Na
    ID = asap(a).name;
    ls = length(asap(a).samecells);
    
    for s = 1:ls                                                            % For each group of registered cells
        figure('Name',[ID,' group ',num2str(s)]);
        reg = asap(a).samecells{s};
        lr = size(reg,1);                                                   % # days a cell was followed
        for r = 1:lr                                                        % For each ROI (day) within this group
            % FIND CORRECT FILE TO LOAD
            dd = reg(r,1);                                                  % Index of imaging day
            day = asap(a).sessions{dd,1};                                   % Actual date
            sess = reg(r,2);                                                % Actual session index
            roi = reg(r,3);                                                 % ROI index
            sessfov = find(asap(a).sessions{dd,3} == sess);                 % Index of that session in the sessions-to-keep
            fov = asap(a).sessions{dd,6}(sessfov);                          % FOV index
            
            % LOAD CORRECT FILE (CANNOT REPLACE WITH POOLED DATA DUE TO SESSION NUMBERING IN DAYS_PV/SOM)
            [asapfile,~] = get_ASAPfile(ID,day,sess);
            [path,videoname] = fileparts(asapfile);
            Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
            modfile = fullfile(path,[videoname,'_Mcells.mat']);
            load(modfile,'Mcells','trials','bins','timepoints','onbins');
            load(Datafile,'V','B');
            disp(asapfile);
            
            % GET RATES AND CELL GROUP
            R = V.R;
            R = squeeze(R(roi,:,:))';                                       % Keep selected roi rates [trials x bins]
            R = (R - mean(R(:,onbins),2)) ./ std(R(:,onbins),[],2);         % ZSCORE THE RATE OVER MODULATION BINS ONLY!!!
            R(isnan(R) | isinf(R)) = 0;
            
            [Ntr,lt] = size(R);
            
            mcells = Mcells(roi,:);                                         % Load Mcell identity of roi
            celltype = mcells(1);                                           % Keep the cell type
            mbin = mcells(2);                                               % and its field bin
         
            % GET PROGRESS AND LICK TIMES
            progress = B.progress;
            licks = B.licks;
            licks = get_licks(licks,timepoints);
            
            %% PLOT FOVs
            FOV = Tiff(fullfile(path,['AVG_FOV',num2str(fov),'.tif']),'r'); % Read the average FOV file
            FOV = read(FOV);
            
            % CROP FOV
            w = 100;
            [~,m] = max(FOV,[],'all','linear');                             % Center FOV around the pixel with max signal
            [x,y] = ind2sub(size(FOV),m);                                   % Turn pixel number to coordinates
            FOV([1:x-w, x+w:end], :) = [];                                  % Keep FOV centered there and within window
            FOV(:, [1:y-w, y+w:end]) = [];

            subplot(5,lr,r);
            imshow(FOV,[min(FOV(:))*0.8, max(FOV(:))]);
            axis square
            set(gca,'Xtick',{},'Ytick',{})
            title(['Day ',num2str(dd)]);
            
            %% PLOT CELL ACTIVITY ACROSS DAYS
            % ORGANIZE TRIAL SEQUENCE
            combos = progress(:,1);                                         % Keep the trial odor combos
            [~,~,combos] = unique(combos);                                  % Order them based on combo
            [combos,trorder] = sort(combos);                                % Sort them by AA AB BA BB, and keep new trial order
            
            R = R(trorder,:);                                               % Sort rates etc
            trials = trials(trorder);
            licks = licks(trorder);
            
            % MEAN RATES BY FIRSTODOR
            Mt = zeros(2,lt);
            St = zeros(2,lt);
            for t = 1:2                                                     % For each type of trial
                Mt(t,:) = mean(R(trials == t,:),1);                         % Keep the cell's mean firing rate
                St(t,:) = SEM(R(trials == t,:));                            % Keep the cell's mean firing rate
            end
            
            % GET SIGNIFICANCE OF MEAN RATES BETWEEN THE TWO ODORS
            onbins = find(onbins);
            p = ones(1,length(onbins));
            for t = 1:length(onbins)
                p(t) = ranksum(squeeze(R(trials == 1,onbins(t))),squeeze(R(trials == 2,onbins(t))));
            end
            [~,p] = fdr(p);
            
            % PLOT EACH TRIAL TYPE
            subplot(5,lr,([1 2 3]*lr +r));
            plot_seq_rates(R,[],bins,timepoints,2);                         % Plot trials split by type (color)
            for tt = 1:3
                k = sum(combos <= tt);
                line([bins(1) timepoints(7)],[1 1]*k+0.5,'color','w','linewidth',2); % Add lines between trial types
            end
            xlim([0 timepoints(5)]);
            for tr = 1:Ntr
                plot(licks{tr},tr*ones(size(licks{tr})),'.w');
            end
            
            if celltype == 1
                line([1 1]*mbin,[1 sum(trials == 1)],'color','r','linestyle','--','linewidth',2); % Add dashed line on its field bin
            elseif celltype == 2
                line([1 1]*mbin,[sum(trials == 1)+1 Ntr],'color','r','linestyle','--','linewidth',2); % Add dashed line on its field bin
            elseif celltype == 0
                line([1 1]*mbin,[1 Ntr],'color','r','linestyle','--','linewidth',2); % Add dashed line on its field bin
            end
            
            subplot(5,lr,4*lr + r); hold on;
            plot_mrate_traces(Mt,St,bins,timepoints,C);                     % Plot odor-average rates
            plot(bins(onbins(p < 0.05)), ones(sum(p < 0.05),1)*(max(Mt(:))+0.2),'.k'); % and significant differences
        end
        count = count + 1;
        
        disp(' ');
    end
end
