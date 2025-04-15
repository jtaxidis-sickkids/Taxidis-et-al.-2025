function Plot_RateCorr_vs_Days(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV
elseif strcmp(INclass,'SOM')
    Days_SOM
end

Na = length(asap);                                                          % total number of animals
% Nr = length([asap.samecells]');                                             % total number of cells followed across days

%% POOL RATES OVER TWO CONSECUTIVE IMAGING SESSIONS OR DO PRE-POST TRAINING COMPARISON
Corr = cell(Na,3);
Ddist = cell(Na,3);

for a = 1:Na                                                                % For each animal
    ID = asap(a).name;
    ls = length(asap(a).samecells);
    
    for s = 1:ls                                                            % For each group of registered cells
        reg = asap(a).samecells{s};                                         % Keep the info of all registered cells
        days = reg(:,1);                                                    % Keep the day indexes
        
        Days = cell(1,2);
        Days{1} = find(days > 2);                                           % Turn to row vector
        last_pre = find(days <= 2,1,'last');                                % Find the last naive imaging day of the cell (keep index in list of days)
        first_post = find(days > 2,1,'first');                              % and the first trained day
        if ~isempty(last_pre) && ~isempty(first_post)                       % If the cell was imaged before and after training
            Days{2} = [last_pre; first_post];                               % Keep that last naive day
        end

        for i = 1:2
            for d1 = 1 : length(Days{i})-1                                  % For each selected registered day (except last)
                % FIND CORRECT FILE TO LOAD
                D1 = Days{i}(d1);
                Sess1 = reg(D1,1);                                          % Actual imaging session number
                day = asap(a).sessions{Sess1,1};                            % Actual date
                sess = reg(D1,2);                                           % Actual video index
                roi = reg(D1,3);                                            % ROI index
                
                % LOAD CORRECT FILE (CANNOT REPLACE WITH POOLED DATA DUE TO SESSION NUMBERING IN DAYS_PV/SOM)
                [asapfile,~] = get_ASAPfile(ID,day,sess);
                [path,videoname] = fileparts(asapfile);
                Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
                modfile = fullfile(path,[videoname,'_Mcells.mat']);
                load(modfile,'onbins');
                load(Datafile,'V');
                disp(asapfile);
                
                R = V.R;
                R = squeeze(R(roi,:,:))';                                   % Keep roi rates [trials x bins]
                R = R(:,onbins);
                R = zscore(R,[],2);                                         % ZSCORE THE RATE OVER MODULATION BINS ONLY!!!
                R(isnan(R) | isinf(R)) = 0;
                R1 = R;
                
                for d2 = d1+1 : length(Days{i})                             % For that day and the next on the session list
                    % FIND CORRECT FILE TO LOAD
                    D2 = Days{i}(d2);
                    Sess2 = reg(D2,1);                                      % Actual imaging session number
                    day = asap(a).sessions{Sess2,1};                        % Actual date
                    sess = reg(D2,2);                                       % Actual video index
                    roi = reg(D2,3);                                        % ROI index
                    
                    % LOAD CORRECT FILE (CANNOT REPLACE WITH POOLED DATA DUE TO SESSION NUMBERING IN DAYS_PV/SOM)
                    [asapfile,~] = get_ASAPfile(ID,day,sess);
                    [path,videoname] = fileparts(asapfile);
                    Datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
                    modfile = fullfile(path,[videoname,'_Mcells.mat']);
                    load(modfile,'onbins');
                    load(Datafile,'V');
                    disp(asapfile);
                    
                    R = V.R;
                    R = squeeze(R(roi,:,:))';                               % Keep roi rates [trials x bins]
                    R = R(:,onbins);
                    R = zscore(R,[],2);                                     % ZSCORE THE RATE OVER MODULATION BINS ONLY!!!
                    R(isnan(R) | isinf(R)) = 0;
                    R2 = R;
                    
                    C = corr(R1',R2','type','Pearson');                     % Correlations between all combinations of trials from the two days
                    C = nanmean(C,'all');                                   % Keep the average across trial combinations
                    
                    Corr{a,i} = cat(1,Corr{a,i},C);                         % Stack the correlation
                    Ddist{a,i} = cat(1,Ddist{a,i}, Sess2-Sess1);            % and the distance between imaging sessions
                end
            end
            disp(' ');
        end
    end
end

%% SET PRE-POST DAY DISTANCE TO 6
for a = 1:Na
    Ddist{a,2} = 6*ones(size(Ddist{a,2}));
end

%% LEAST SQUARES ESTIMATE FOR DISTRIBUTION
C1 = cell2mat(Corr(:,1));
D1 = cell2mat(Ddist(:,1));
C2 = cell2mat(Corr(:,2));
D2 = cell2mat(Ddist(:,2));

coefficients = polyfit(D1, C1, 1);                    % Equivalent to: [ones(size(D)), D]\C
Cfitted = polyval(coefficients, [0; D1; max(D1)+1]);
pfit = coefTest(fitlm(D1,C1,'linear'));

%% SPEARMAN CORRELATION OF DISTRIBUTIONS ACROSS DAY-DISTANCE
[R,P] = corr(D1,C1,'type','Spearman')

%% PLOT POOLED RATES OVER TWO CONSECUTIVE DAYS
figure; 
subplot(1,3,1:2);hold on;
plot(D1 + randn(size(D1))*0.1,C1,'ok','MarkerFaceColor','k')
plot([0; D1; max(D1)+1] , Cfitted,'-k')

plot(D2 + randn(size(D2))*0.1,C2,'ok')

xlabel('Distance between sessions');
ylabel('Firing rates average correlation');
title(num2str(pfit));

%% COMPARE TRAINED VS PRE-POST
subplot(133);hold on;
p = ranksum(C1,C2);

cc = table(C1',C2','VariableNames',{'trained','pre-post'});
violinplot(cc);
plot_significance(p,1,2,C1,C2)
title(num2str(p));
ylabel('Correlations');





