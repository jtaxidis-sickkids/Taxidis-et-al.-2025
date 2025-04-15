function Compare_Behavior(sessions)

Poff = nan(length(sessions),4);
Pon = Poff;

for s = 1:length(sessions)
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'opto','progress');
    
    outcomes = progress(:,2);                                               % Get the outcome of each trial
    correct = outcomes == 1 | outcomes == 4;                                % Find corrects
    outcomes(correct) = 1;                                                  % Set those as = 1
    outcomes(~correct) = 0;                                                 % And errors as = 0
    
    for i = 1:3                                                             % For each opto protocol with intermixed light on/off trials
        koff = opto(:,1) == i & opto(:,2) == 0;                             % Keep all the light-off trials of that protocol
        kon = opto(:,1) == i & opto(:,2) == 1;                              % and the light-on
        
        Poff(s,i) = sum(outcomes(koff)) / sum(koff);                        % Compute corresponding performances
        Pon(s,i) = sum(outcomes(kon)) / sum(kon);
    end
    
    % REPEAT FOR NO-OPTO vs FULL-ODOR OPTO
    koff = opto(:,1) == 0;                                                  % Repeat by comparing trials with no opto
    kon = opto(:,1) == 4;                                                   % with full-odor opto
    Poff(s,4) = sum(outcomes(koff)) / sum(koff);
    Pon(s,4) = sum(outcomes(kon)) / sum(kon);
end

Poff = Poff*100;                                                            % Turn to % percentage
Pon = Pon*100;

figure;
for i = 1:4
    subplot(2,4,i); hold on;
    plot([1 2],[Poff(:,i), Pon(:,i)],'o-b');
    errorbar([0.9 2.1],[nanmean(Poff(:,i)) nanmean(Pon(:,i))],[SEM(Poff(:,i)) SEM(Pon(:,i))],'ok-')
    set(gca,'Xtick',[1,2],'XTickLabel',{'no opto';'opto'});
    if i == 1, ylabel('Performance per day'); end
    
    [~,pP] = ttest(Poff(:,i),Pon(:,i))
    plot_significance(pP,0.9,2.1,nanmean(Poff(:,i)), nanmean(Pon(:,i)))
end


%% REPEAT BUT PER TRIAL BLOCK
Poff = cell(1,4);
Pon = Poff;
counter = [1 1 1 1];

for s = 1:length(sessions)
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'opto','progress');
    
    outcomes = progress(:,2);                                               % Get the outcome of each trial
    correct = outcomes == 1 | outcomes == 4;                                % Find corrects
    outcomes(correct) = 1;                                                  % Set those as = 1
    outcomes(~correct) = 0;                                                 % And errors as = 0
    
    Ntr = length(outcomes);
    
    for i = 1:3                                                             % For each opto protocol with intermixed light on/off trials
        Nblocks = sum(opto(:,1) == i)/20;
        block = repmat([1:Nblocks],20,1);
        block = block(:);
        Block = zeros(1,Ntr)';
        Block(opto(:,1) == i) = block;
        
        for j = 1:Nblocks
            koff = opto(:,1) == i & opto(:,2) == 0 & Block == j;                             % Keep all the light-off trials of that protocol
            kon  = opto(:,1) == i & opto(:,2) == 1 & Block == j;                              % and the light-on
            
            Poff{i}(counter(i)) = sum(outcomes(koff)) / sum(koff);                        % Compute corresponding performances
            Pon{i}(counter(i)) = sum(outcomes(kon)) / sum(kon);
            
            counter(i) = counter(i) + 1;
        end       
    end
    
    % REPEAT FOR NO-OPTO vs FULL-ODOR OPTO   
    koff = opto(:,1) == 0;                                                  % Repeat by comparing trials with no opto
    kon = opto(:,1) == 4;                                                   % with full-odor opto
    Poff{4}(s) = sum(outcomes(koff)) / sum(koff);    
    Pon{4}(s) = sum(outcomes(kon)) / sum(kon);
end

for i = 1:4
    Poff{i} = Poff{i}'*100;                                                 % Turn to columns and to percentage
    Pon{i} = Pon{i}'*100;
end

for i = 1:4
    subplot(2,4,4+i); hold on;
    plot([1 2],[Poff{i}, Pon{i}],'o-b');
    errorbar([0.9 2.1],[nanmean(Poff{i}) nanmean(Pon{i})],[SEM(Poff{i}) SEM(Pon{i})],'ok-')
    set(gca,'Xtick',[1,2],'XTickLabel',{'no opto';'opto'});
    if i == 1, ylabel('Performance per block'); end
    
    [~,pP] = ttest(Poff{i},Pon{i})
    plot_significance(pP,0.9,2.1,nanmean(Poff{i}), nanmean(Pon{i}))
end


