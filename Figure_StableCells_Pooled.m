function Figure_StableCells_Pooled(INclass,plot_flag)

if strcmp(INclass,'PV')
    Days_PV
    maxd = 7;
elseif strcmp(INclass,'SOM')
    Days_SOM
    maxd = 8;
end

Na = length(asap);                                                          % total number of animals
% Nr = length([asap.samecells]');                                             % total number of cells followed across days

%% POOL RATES OVER TWO CONSECUTIVE IMAGING SESSIONS OR DO PRE-POST TRAINING COMPARISON
MR = [];
F = [];
ISMCELL = [];

for a = 1:Na                                                                % For each animal
    ID = asap(a).name;
    ls = length(asap(a).samecells);
    
    for s = 1:ls                                                            % For each group of registered cells
        reg = asap(a).samecells{s};                                         % Keep the info of all registered cells 
        days = reg(:,1);                                                    % Keep the day indexes
        
        selected_days = [];
        if strcmp(plot_flag,'all')
            selected_days = find(days);                                     % Keep all imaging days (their indexes in the list of days)
            
        elseif strcmp(plot_flag,'trained')
            selected_days = find(days > 2);                                 % Keep only trained imaging days (their indexes in the list of days)
            
        elseif strcmp(plot_flag,'pre_post')
            last_pre = find(days <= 2,1,'last');                            % Find the last naive imaging day of the cell (keep index in list of days)
            first_post = find(days > 2,1,'first');                          % and the first trained day
            if ~isempty(last_pre) && ~isempty(first_post)                   % If the cell was imaged before and after training
                selected_days = last_pre;                                   % Keep that last naive day                
                selected_days = [selected_days; nan];                       % (Add a NaN for the for-loop below to run)
            end
        end
        selected_days = selected_days';                                     % Turn to row vector
        
        for d = selected_days(1:end-1)                                    % For each selected registered day (except last)        
            mr = [];
            mbin = [];
            ismcell = [];

            for dd = [d d+1]                                                % For that day and the next on the session list
                % FIND CORRECT FILE TO LOAD
                sessionindx = reg(dd,1);                                    % Actual imaging session number
                day = asap(a).sessions{sessionindx,1};                      % Actual date
                sess = reg(dd,2);                                           % Actual video index
                roi = reg(dd,3);                                            % ROI index
                
                % LOAD CORRECT FILE (CANNOT REPLACE WITH POOLED DATA DUE TO SESSION NUMBERING IN DAYS_PV/SOM)
                [asapfile,~] = get_ASAPfile(ID,day,sess);
                [path,videoname] = fileparts(asapfile);
                Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
                modfile = fullfile(path,[videoname,'_Mcells.mat']);
                load(modfile,'Mcells','bins','onbins','timepoints');
                load(Datafile,'V');
                disp(asapfile);
                
                R = V.R;
                R = squeeze(R(roi,:,:))';                                   % Keep roi rates [trials x bins]
                R = (R - mean(R(:,onbins),2)) ./ std(R(:,onbins),[],2);     % ZSCORE THE RATE OVER MODULATION BINS ONLY!!!
                R(isnan(R) | isinf(R)) = 0;
                
                mr = cat(3,mr,mean(R,1));                                   % Stack the two average rates in the 3rd dim
                mbin = [mbin, Mcells(roi,2)];                               % Stack the two maxbins next to each other
                ismcell = [ismcell, ~isnan(Mcells(roi,1))];                 % Same for flag showing if cell has any field
            end
            
            MR = [MR; mr];                                                  % Store pooled data from all cells and all pairs of days
            F = [F; mbin];
            ISMCELL = [ISMCELL; ismcell];
            
            disp(' ');
        end
    end
end

%% SORT POOLED DATA BASED ON DAY 1
[~,order] = sort(F(:,1));
F = F(order,:);
MR = MR(order,:,:);
ISMCELL = ISMCELL(order,:);
ISMCELL = logical(ISMCELL);

%% PLOT POOLED RATES OVER TWO CONSECUTIVE DAYS
figure;
for d = 1:2
    subplot(1,2,d);
    plot_seq_rates(MR(:,:,d),[],bins,timepoints,1);       % Plot pooled sequences
    plot(F(:,d),1:size(F,1),'.k')
    plot(F(ISMCELL(:,d),d),find(ISMCELL(:,d)),'oy','MarkerFaceColor','y')
    xlim([0, timepoints(5)]);
    title(['Day ',num2str(d)]);
end

