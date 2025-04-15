function [MRp,Fp] = Pool_sequences(sessions,optnum,optoON)

celltypes = {'PY';'IN'};
TP = [1 2 7 8 9 11 11];                                     

%% LOAD ALL PROCESSED SEQUENCES IN THE CHOSEN FILES
for s = 1:length(sessions)                                           
    dirname = ['../ProcessedData/',sessions{s}];
    DL(s) = load(fullfile(dirname,'CA1data.mat'),'R','opto','progress','bins');
    DPI(s) = load(fullfile(dirname,'PY_IN.mat'),'PY_IN');
    DM(s) = load(fullfile(dirname,'Mcells.mat'),'Mcells','trials');
end

ctypes = 2;
ttypes = 2;

ls = length(sessions);

bins = DL(1).bins;                                        
bins(end) = [];
lb = length(bins);

%% POOL AND SORT M-CELLS
MRp = cell(ctypes,ttypes,ttypes);
Fp = cell(ctypes,ttypes);
% SIp = cell(ctypes,ttypes);

for ct = 1:ctypes                                                           % For each cell type
    for tt1 = 1:ttypes                                                      % For each trial type
        Mc = cell(ls,1);
        F = cell(ls,1);
        P = cell(ls,1);
%         SI = cell(ls,1);
        for s = 1:ls                                                        % For each sesssion     
            mcells = DM(s).Mcells{ct,tt1};
            if isempty(mcells); mcells = zeros(0,4); end
            
            Mc{s} = mcells(:,1);                                            % Keep the mcells
            F{s} = mcells(:,2);                                             % and their fields
            P{s} = mcells(:,3);
%             SI{s} = mcells(:,4);
        end
        
        F = cell2mat(F);                                                    % Concatenate fields
        [F,order] = sort(F);                                                % Sort fields
        Fp{ct,tt1} = F;                                                     % and save
        
%         SI = cell2mat(SI);                                                    % Concatenate fields
%         SI = SI(order);
%         SIp{ct,tt1} = SI;
        
        for tt2 = 1:ttypes
            MR = cell(ls,1);
            for s = 1:ls                                                    % For each session
                trials = DL(s).progress;                                    % Keep all trials
                trials = floor(trials(:,1)/10);                             % Get their first odor
                opto = DL(s).opto;                                          % and the opto protocols
                
                R = DL(s).R;                                                % Take all firing rates    
                R = R(DPI(s).PY_IN == ct,:,:);                              % KEEP ONLY CELLS OF THAT CELLTYPE
                R = R(Mc{s},:,:);                                           % Keep only mcells
                         
                if optnum == 0
                    optotrials = (opto(:,2) == optnum & trials == tt2);      % KEEP ALL TRIALS WITHOUT LIGHT AND THAT ODOR
                else
                    optotrials = (opto(:,1) == optnum & opto(:,2) == optoON & trials == tt2);      % KEEP ONLY TRIALS OF THAT OPTO AND THAT ODOR
                end
                R = R(:,:,optotrials);
                
                R = smoothdata(R,2,'gaussian',30);                           % SMOOTH EACH TRIAL (SAME AS IN MCELL CALCULATION)
                
                MR{s} = nanmean(R,3);                                          % Keep mean signals of all cells of that type over trials of that type               
%                 MR{s} = MR{s} ./ repmat(P{s},1,lb); % mean(MR{s}(:,bins < 1),2);%             % Normalize by field signal
                MR{s} = MR{s} ./ max(MR{s}(:,bins>1 & bins <7),[],2);%             % Normalize by max odor/delay rate
%                 MR{s} = MR{s} ./ mean(MR{s}(:,bins<1),2);%             % Normalize by mean baseline
            end
            MR = cell2mat(MR);                                              % Concatenate rates of all cells over all days
            MR = MR(order,:);                                               % Sort by field time           
            MRp{ct,tt1,tt2} = MR;                                           % Store
        end
    end
end

%% PLOT AVERAGE SIGNALS (COLOR OR LINES)
for ct = 1:ctypes                                                       % For each cell type
    figure('Name',['Pooled ',celltypes{ct},' M-cells']);
    for tt1 = 1:ttypes                                                  % For each trial type
        for tt2 = 1:ttypes                                              % For each trial type (again)
            subplot(ttypes,ttypes,(tt1-1)*ttypes + tt2); hold on;
            plot_seq_rates(MRp{ct,tt1,tt2},[],bins,TP,1); % Plot pooled sequences
        end
        subplot(ttypes,ttypes,(tt1-1)*ttypes + 1);
        ylabel(['Sequence-',num2str(tt1)]);
    end
end


%% PLOT ALL MCELLS ON THEIR PREFERRED VS NON-PREFERRED TRIALS
for ct = 1:ctypes                                                           % For each cell type
    figure('Name',['Pooled ',celltypes{ct},' M-cells, Preferred trials']);
    F = [Fp{ct,1}; Fp{ct,2}];
    [F,order] = sort(F);
    
    Mp = [MRp{ct,1,1}; MRp{ct,2,2}];
    Mp = Mp(order,:);
    
    subplot(121);
    plot_seq_rates(Mp,[],bins,TP,1);       % Plot pooled sequences
    title('Preferred trials');
    xlim([0 9]);
    
    plot(F,1:length(F),'ow','MarkerFaceColor','w','Markersize',3);
    
    z = find(F >= TP(2),1,'first');
    if ~isempty(z)
        line([bins(1) bins(lb)], z*[1 1],'Color','w','LineStyle','--','LineWidth',2); % Add line where time-cells start
    end
    
    Mnp = [MRp{ct,1,2}; MRp{ct,2,1}];
    Mnp = Mnp(order,:);
    
    subplot(122);
    plot_seq_rates(Mnp,[],bins,TP,1);      % Plot pooled sequences
    title('Non-preferred trials');
    xlim([0 9]);

    z = find(F >= TP(2),1,'first');
    if ~isempty(z)
        line([bins(1) bins(end)], z*[1 1],'Color','w','LineStyle','--','LineWidth',2);
    end
end




