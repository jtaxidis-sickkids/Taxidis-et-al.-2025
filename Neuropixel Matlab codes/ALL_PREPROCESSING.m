
clear
sessions = {... 
            % NAIVE
            'ZD028_0822', 'WT';
            'ZD028_0824', 'WT';
%             'ZD029_0823', 'WT';  % Clustered spearately 
%             'ZD029_0824', 'WT';  % Clustered spearately 
            'CD234_0825', 'PV';
            'CD234_0826', 'PV';
            'CD235_0826', 'PV';
            'ZD031_0827', 'SST';
            'ZD031_0828', 'SST'; 
            'ZD032_0827', 'SST';
            'ZD032_0828', 'SST';
            
%             TRAINED
            'ZD028_0904', 'WT';
%             'ZD029_0904', 'WT'; % Clustered spearately  
            'CD234_0905', 'PV';
            'CD234_0907', 'PV';
            'CD235_0905', 'PV';
            'CD235_0906', 'PV';
            'ZD031_0906', 'SST';
            'ZD031_0907', 'SST';
            'ZD031_0908', 'SST';
};%

ls = size(sessions,1);

session = 'CD234_0905';  %

% SET OPTO PROTOCOLS
[sessiontypes, Nprobes] = Set_Sessiontypes(session)

%% LOAD NEUROPIXEL/WINEDR/MATLAB DATA, KEEP GOOD UNITS, AND SAVE IN ONE FILE
Preprocess(session,sessiontypes,Nprobes)

%% FIND THE CA1-LOCATED CHANNEL
CA1ch = Find_CA1(session);

%% KEEP CA1 DATA ONLY AND ORGANIZE DATA INTO TRIALS
Organize_trials(session,80)

%% ANALYZE WAVEFORMS
Waveform_Analysis(session)

%% SPLIT INTO PY-IN
Unit_Clustering(sessions(:,1))

%% FIND SEQUENCE CELLS
for s = 1:ls
    Get_Modulation(sessions{s,1});
end
