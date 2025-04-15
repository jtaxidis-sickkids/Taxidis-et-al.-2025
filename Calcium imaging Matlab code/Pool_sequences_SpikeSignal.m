function [MP,MNP,FF] = Pool_sequences_SpikeSignal(days,crit)

celltypes = {'PY';'IN'};

%% LOAD ALL PROCESSED SEQUENCES IN THE CHOSEN FILES
files = [];
for dd = 1:length(days)                                                     % For each deay
    files = [files; subdir(fullfile(days{dd},['Modulation*',crit,'_ASAP.mat']))]; % Get all Modulation files (all delays)
end

for d = 1:length(files)                                                     % For each Modulation file
    parent_folder = fileparts(files(d).name);
    [parent_folder,day] = fileparts(parent_folder);
    disp([parent_folder,'    ',day]);
    
    DL(d) = load(files(d).name);                                            % Load it
    delays(d) = DL(d).timepoints(3) - DL(d).timepoints(2);                  % and COMPUTE the delay
end
[delays,order] = sort(delays);                                              % Sort the delays of ALL files (not unique)
DL = DL(order);                                                             % and rearrange all files

uniqdel = unique(delays);                                                   % Get the different delays
ldel = length(uniqdel);

[ctypes,ttypes] = size(DL(1).MR);                                           % Get number of cell types and trial types
t = DL(1).time;                                                             % and keep all info needed
t1 = t(1);
dt = t(2)-t(1);

%% POOL AND SORT M-CELLS
MSp = cell(ctypes,ttypes,ttypes,ldel);
SSp = cell(ctypes,ttypes,ttypes,ldel);
Fp = cell(ctypes,ttypes,ldel);
TPp = [];
for d = 1:ldel                                                              % For each delay
    dsess = find(delays == uniqdel(d));                                     % Find all files of that delay
    ls = length(dsess);
    TPp(d,:) = DL(dsess(1)).timepoints;                                     % keep their timepoints
    ontime{d} = DL(dsess(1)).ontime;                                        % and modulation bins

    for ct = 1:ctypes                                                       % For each cell type
        for tt1 = 1:ttypes                                                  % For each trial type
            Mc = cell(ls,1);
            F = cell(ls,1);
            P = cell(ls,1);
            lt = zeros(ls,1);
            for s = 1:ls                                                    % For each file with that delay
                if isempty(DL(dsess(s)).Mcells{ct,tt1}); DL(dsess(s)).Mcells{ct,tt1} = zeros(0,3); end
          
                Mc{s} = DL(dsess(s)).Mcells{ct,tt1}(:,1);                   % Keep the mcells
                F{s} = DL(dsess(s)).Mcells{ct,tt1}(:,2);                    % and their fields
                P{s} = DL(dsess(s)).Mcells{ct,tt1}(:,3);
                lt(s) = length(DL(dsess(s)).time);                          % and the number of bins (changes between days depending on shortest trial of that day)
            end
            minlt(d) = min(lt);                                             % Store the minimum number of bins over all days
            
            F = cell2mat(F);                                                % Concatenate fields
            [F,order] = sort(F);                                            % Sort fields
            Fp{ct,tt1,d} = F;                                               % and save

            for tt2 = 1:ttypes
                MS = cell(ls,1);
                SS = cell(ls,1);
                for s = 1:ls                                                % For each file with that delay
                    if size(DL(dsess(s)).MR,1) == 1                         % For mice with only PY cells (ctypes ==1)
                        DL(dsess(s)).MR = [DL(dsess(s)).MR; {zeros(0,minlt(d)),zeros(0,minlt(d))}];
                        DL(dsess(s)).SR = [DL(dsess(s)).SR; {zeros(0,minlt(d)),zeros(0,minlt(d))}];
                    end
                    
                    MS{s} = DL(dsess(s)).MR{ct,tt2};                        % Keep corresponding mean signals of all cells of that type over trials of that type
                    SS{s} = DL(dsess(s)).SR{ct,tt2};                        % Keep corresponding SD signals
                    
                    MS{s} = MS{s}(Mc{s},:);                                 % Keep only mcells signals
                    SS{s} = SS{s}(Mc{s},:);                                 % Same for SDs
                    
                    MS{s} = MS{s}(:,1:minlt(d));                            % Keep only the signal over the min number of bins
                    SS{s} = SS{s}(:,1:minlt(d));                            % Same for SDs
                    
                    MS{s} = MS{s} ./ repmat(P{s},1,minlt(d));               % Normalize by field signal
                    SS{s} = SS{s} ./ repmat(P{s},1,minlt(d));
                end            
                MS = cell2mat(MS);                                          % Concatenate rates of all cells over all days
                SS = cell2mat(SS);                                          % and SDs
                
                MS = MS(order,:);                                           % Sort by field time
                SS = SS(order,:);
             
                MSp{ct,tt1,tt2,d} = MS;                                     % Store 
                SSp{ct,tt1,tt2,d} = SS;                                           
            end
        end
    end
end

%% PLOT AVERAGE SIGNALS (COLOR OR LINES)
for d = 1:ldel                                                              % For each delay
    for ct = 1:ctypes                                                       % For each cell type
        figure('Name',['Pooled ',celltypes{ct},' M-cells, Delay = ',num2str(uniqdel(d))]);
        for tt1 = 1:ttypes                                                  % For each trial type
            for tt2 = 1:ttypes                                              % For each trial type (again)
                subplot(ttypes,ttypes,(tt1-1)*ttypes + tt2); hold on;
                plot_seq_rates(MSp{ct,tt1,tt2,d},SSp{ct,tt1,tt2,d},[],ontime{d},t(1:minlt(d)),TPp(d,:),2); % Plot pooled sequences
                title(['Trials: ',crit,'-',num2str(tt2)]);
            end
            subplot(ttypes,ttypes,(tt1-1)*ttypes + 1);
            ylabel(['Sequence-',num2str(tt1)]);
        end
    end
end

%% PLOT ALL MCELLS ON THEIR PREFERRED VS NON-PREFERRED TRIALS
for d = 1:ldel                                                              % For each delay   
    for ct = 1:ctypes                                                       % For each cell type
        figure('Name',['Pooled ',celltypes{ct},' M-cells, Preferred trials. Delay = ',num2str(uniqdel(d))]);
        F = [Fp{ct,1,d}; Fp{ct,2,d}];
        [F,order] = sort(F);
        
        Mp = [MSp{ct,1,1,d}; MSp{ct,2,2,d}];
        Mp = Mp(order,:);
        
        subplot(121);
        plot_seq_rates(Mp,[],[],ontime,t1:dt:dt*minlt(d),TPp(d,:),2);       % Plot pooled sequences
        title('Preferred trials');
       
        z = find(F >= TPp(d,2),1,'first');
        line([t1 dt*minlt(d)], z*[1 1],'Color','w','LineStyle','--','LineWidth',2); % Add line where time-cells start
        
        Mnp = [MSp{ct,1,2,d}; MSp{ct,2,1,d}];
        Mnp = Mnp(order,:);
        
        subplot(122);
        plot_seq_rates(Mnp,[],[],ontime,t1:dt:dt*minlt(d),TPp(d,:),2);      % Plot pooled sequences
        title('Non-preferred trials');
       
        z = find(F >= TPp(d,2),1,'first');
        line([t1 dt*minlt(d)], z*[1 1],'Color','w','LineStyle','--','LineWidth',2);

        if uniqdel(d) == 5 && ct == 2           % !!!!!!!!!!!!
          MP = Mp;
          MNP = Mnp;
          FF = F;
        end
    end
end



