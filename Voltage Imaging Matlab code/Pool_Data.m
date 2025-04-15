function Pool_Data(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
    maxd = 7;
elseif strcmp(INclass,'SOM')
    Days_SOM;
    maxd = 8;
end

Na = length(asap);

DFF = cell(Na,maxd);
R = cell(Na,maxd);
SP = cell(Na,maxd);
MC = cell(Na,maxd);
NMC = cell(Na,maxd);
TR = cell(Na,maxd);
AMP = cell(Na,maxd);
PH = cell(Na,maxd);
MDIR = cell(Na,maxd);
LVEC = cell(Na,maxd);
MOT = cell(Na,maxd);
PROG = cell(Na,maxd);

Fs = 1000;
time = 0 : 1/Fs : (11*Fs-1)/Fs;

%% LOAD ALL PROCESSED SEQUENCES IN THE CHOSEN FILES
for a = 1:length(asap)
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1:length(days)                                               % For each selected day
        MC{a,d} = {};
        NMC{a,d} = {};
        TR{a,d} = {};
        DFF{a,d} = {};
        R{a,d} = {};
        SP{a,d} = {};
        AMP{a,d} = {};
        PH{a,d} = {};
        MDIR{a,d} = {};
        LVEC{a,d} = {};
        MOT{a,d} = {};
        PROG{a,d} = {};
        
        for s = 1:length(ses{d})                                    % for each selected session that day
            [asapfile,~] = get_ASAPfile(ID,days{d},ses{d}(s));
            [path,videoname] = fileparts(asapfile);
            disp(asapfile);
            modfile = fullfile(path,[videoname,'_Mcells.mat']);
            load(modfile,'Mcells','trials','bins','timepoints','onbins');
            
            MC{a,d} = cat(1,MC{a,d},Mcells);
            TR{a,d} = cat(1,TR{a,d},trials);

            modfile = fullfile(path,[videoname,'_Neg_Mcells.mat']);
            load(modfile,'Mcells');
            
            NMC{a,d} = cat(1,NMC{a,d},Mcells);         
            
            datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
            load(datafile,'V','F','M','B');           
            
            DFF{a,d} = cat(1,DFF{a,d},V.DespDFF);
            R{a,d} = cat(1,R{a,d},V.R);
            SP{a,d} = cat(1,SP{a,d},{V.SP});
            
            AMP{a,d} = cat(1,AMP{a,d},F.Ampl);
            PH{a,d} = cat(1,PH{a,d},F.Phase);
            MDIR{a,d} = cat(1,MDIR{a,d},F.mDir);
            LVEC{a,d} = cat(1,LVEC{a,d},F.lVec);
            
            MOT{a,d} = cat(1,MOT{a,d},M.MOT); 
            
            PROG{a,d} = cat(1,PROG{a,d},B.progress); 
        end
    end
end

%% SAVE
save(['PooledData_',INclass,'.mat'], 'MC','NMC','TR','DFF','R','SP','AMP','PH','MDIR','LVEC','MOT','PROG',...
                                     'time','timepoints','bins','onbins')
