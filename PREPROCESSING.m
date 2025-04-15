
%% VOLPY PROCESSING OF ALL SESSIONS
Days_PV;

for a = 1:length(asap)
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    for d = 1:length(days)
        [~,asapfiles] = get_ASAPfile(ID,days{d},1);
        for s = 1:length(asapfiles)
            VOLPY_processing(ID,days{d},s)
        end
    end
end

%% PLOT VOLPY PROCESSING OF SELECTED SESSIONS
Days_SOM;

for a = 4%:length(asap)
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 3%:length(days)
        for s = 1:length(ses{d})
            Plot_VOLPY_processing(ID,days{d},ses{d}(s))
        end
    end
end

%% PROCESS WAVEFORMS, BURSTINESS, SPIKE PHASE LOCKING AND MOTION
Days_PV;

for a = 1%1:length(asap)
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1%1:length(days)     
        for s = 1:length(ses{d})
            Unit_analysis(ID,days{d},ses{d}(s));
%             Plot_UnitAnalysis(ID,days{d},ses{d}(s));
        end
    end
end
disp('Done');

%% FIND MODULATED CELLS
Days_PV;

for a = 1:length(asap)
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1:length(days)
        for s = 1:length(ses{d})
            Get_Modulation(ID,days{d},ses{d}(s));             % MOSTLY "ODOR-CELLS"
        end
    end
end
disp('Done');

%% FIND NEGATIVELY MODULATED CELLS
Days_PV;

for a = 1:length(asap)
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1:length(days)
        for s = 1:length(ses{d})
            Get_Negative_Modulation(ID,days{d},ses{d}(s));      % MOSTLY "TIME-CELLS'
        end
    end
end
disp('Done');

%% POOL DATA ACROSS ALL SESSIONS
Pool_Data('PV');
Pool_Data('SOM');





