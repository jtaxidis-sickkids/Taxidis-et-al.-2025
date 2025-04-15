
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

% session = 'ZD028_0904';  % one good cell losing it s odorB field and good examples of hyperpol. spiking
% session = 'CD234_0905';  % one good odorA odor cell and one timecell 
% session = 'ZD031_0907';  % good cell released by inhibition at light ON & cell inhibited by odors
% session = 'ZD028_0822';  % good cell inhibited by odors
% session = 'ZD032_0827'; % Good opto effects (Same for next day)
% session = 'ZD028_0904'; % good cells inhibited by odor
% session = 'ZD028_0822'; 
% session = 'CD234_0907';

% OPTO EFFECTS ON TIME CELLS (EXAMPLES) 
session = 'ZD031_0827';
session = 'ZD031_0828';
session = 'ZD032_0827'; %best
% session = 'CD234_0905';

% SET OPTO PROTOCOLS
[sessiontypes, Nprobes] = Set_Sessiontypes(session)

%% COMPARE OPTO
 Compare_Opto_Rates(sessions(:,1))           % Analyzes PV and SST-Cre separately

%% PLOT SPIKES, RATES AND FIELDS OF SESSSION
if strcmp(sessiontypes{end-1}, 'no light')                      
    Plot_Rates(session)
%     Plot_Rates_byOdor(session)
else
    Plot_Rates_opto(session)    
%     Plot_Rates_opto_rearranged(session)
end

%% PLOT POOLED RATES FROM PY AND IN
optocell = 'PV';

k = strcmp(sessions(:,2),optocell);
sessionstemp = sessions(k,1);

% sessionstemp = sessions;          % Uncomment to use all sessions

MR = cell(2,2);
[MR(1,:), order] = Plot_Rates_pooled(sessionstemp,0,0,[]);
[MR(2,:), ~]     = Plot_Rates_pooled(sessionstemp,4,1,order);

% COMPARE POOLED RATES WITHvsWITHOUT OPTO
dirname = ['../ProcessedData/',sessions{1,1}];
load(fullfile(dirname,'CA1data.mat'),'bins');            
bins = bins(bins < 11);

Compare_Rates_Pooled(MR,bins);

%% PLOT POOLED SEQUENCES
optocell = 'PV';

k = strcmp(sessions(:,2),optocell);
sessionstemp = sessions(k,1);
% sessionstemp = sessions;          % Uncomment to use all sessions

[MR_no_opto,F] = Pool_sequences(sessionstemp,0,0);  
[MR_opto,   ~] = Pool_sequences(sessionstemp,1,1);
 
% COMPARE SEQUENCE RATES AND SI SHIFTS WITHvsWITHOUT OPTO
dirname = ['../ProcessedData/',sessions{1,1}];
load(fullfile(dirname,'CA1data.mat'),'bins');   
bins = bins(bins < 7);                         

[SI1,SI2] = Compare_Mcell_MR_SI(MR_no_opto,MR_opto,F,bins,'odor');
[SI1,SI2] = Compare_Mcell_MR_SI(MR_no_opto,MR_opto,F,bins,'delay');
[SI1,SI2] = Compare_Mcell_MR_SI(MR_no_opto,MR_opto,F,bins,'all');           % Use only the odor vs delay SI shift 

%% CHECK LFP THETA RESETTING
optocell = 'SST';

k = strcmp(sessions(:,2),optocell);
sessionstemp = sessions(k,1);

LFP_theta_reset(sessionstemp);

%% CHECK LOCOMOTION
Motion_analysis(sessions(11:14,1))

%% COMPARE BEHAVIOR WITH OPTO
sessionstemp = sessions(11:14,1); % USE ALL PV TRAINED SESSIONS AND WITH OPTO
sessionstemp = sessions(15:17,1); % USE ALL SST TRAINED SESSIONS AND WITH OPTO
% sessionstemp = sessions(11:17,1); % USE ALL TRAINED SESSIONS AND WITH OPTO

Compare_Behavior(sessionstemp);

