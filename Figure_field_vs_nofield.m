function Figure_field_vs_nofield(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
elseif strcmp(INclass,'SOM')
    Days_SOM;
end

%% POOL ISI FROM INDIVIDUAL FILES AND COMPUTE MEDIAN ISI
ISI = {};
ROI = [];
for a = 1:length(asap)                                                                % For each selected animal
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1:length(days)                                               % For each selected day
        for s = 1:length(ses{d})                                    % for each selected session that day
            asapfile = get_ASAPfile(ID,days{d},ses{d}(s));
            [path,videoname] = fileparts(asapfile);
            datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
            load(datafile,'S','V');                                             % Load the Theta-Motion analysis file
            disp(asapfile);
          
            for r = 1:size(S.ISI,1)                                                     % For each selected ROI
                v = cell2mat(S.ISI(r,:));
                ISI = cat(1,ISI,v);
                
                roi = squeeze(V.ROI(r,:,:));
                ROI = cat(1,ROI,sum(roi(:)));   % count ROI pixels
            end 
        end
    end
end
mISI = cellfun(@median, ISI);               % GET EACH CELL's MEDIAN ISI
mISI = mISI * 1000;                         % Turn to msec

%% LOAD POOLED DATA AND STACK
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'R','AMP','MDIR','LVEC','MOT','MC');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    R{1,5} = [];     
    R{2,3} = [];     
    
    AMP{1,5} = [];      
    AMP{2,3} = [];
    
    MDIR{1,5} = [];      
    MDIR{2,3} = [];
        
    LVEC{1,5} = [];      
    LVEC{2,3} = [];
    
    MOT{1,5} = [];     
    MOT{2,3} = [];     
    
    MC{1,5} = [];      
    MC{2,3} = []; 
end

R = vertcat(R{:});
AMP = vertcat(AMP{:});
MDIR = vertcat(MDIR{:});
LVEC = vertcat(LVEC{:});
MOT = vertcat(MOT{:});
MC = vertcat(MC{:});

ls = length(R);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
Rall = {};
AMPall = {};
MDIRall = [];
LVECall = [];
MOTall = {};
MCall = [];

for i = 1:ls
    for r = 1:size(R{i},1)
        Rall{count,1} = squeeze(R{i}(r,:,:))';                          % Keep only that cell's Rates [trials x time]
        MOTall{count,1} = MOT{i};                                           % (Same dimensions)
        AMPall{count,1} = squeeze(AMP{i}(r,:,:,2))';                          % Keep only that cell's AMPLITUDE OVER THETA [trials x time]
        MDIRall(count,1) = MDIR{i}(r,2);                          % Keep only that cell's AMPLITUDE OVER THETA [trials x time]
        LVECall(count,1) = LVEC{i}(r,2);                          % Keep only that cell's AMPLITUDE OVER THETA [trials x time]
        MCall(count,:) = MC{i}(r,:);    
        
        count = count + 1;
    end
end
R = Rall;
AMP = AMPall;
LVEC = LVECall;
MDIR = MDIRall;
MOT = MOTall;
MC = MCall;

clear Rall AMPall MDIRall LVECall MOTall MCall 

ls = count-1;

%% SMOOTH RATES AND MOTION AND DOWNSAMPLE MOTION
% for i = 1:ls
%     MOT{i} = MOT{i}(1:100:end,:);                                           % Downsample to have same vector length as rates
%     for tr = 1:size(R{i},2)
%         R{i}(:,tr) = smooth(R{i}(:,tr),5,'moving');
%     end
%     R{i} = R{i}'; % Turn to [time x trials]
% end

for i = 1:ls
    R{i} = smoothdata(R{i},2,'movmean',5);
    R{i} = R{i}'; % Turn to [time x trials]
    
    MOT{i} = smoothdata(MOT{i},1,'movmean',50);
    MOT{i} = MOT{i}(1:100:end,:);                                           % Downsample to have same vector length as rates
end

%% GET MEAN LOCOMOTION-RATES CORRELATIONS (SPEED SCORE)
Corr = nan(ls,1);
for i = 1:ls
    C = corr(R{i},MOT{i},'type','Pearson');                                 % Correlation of rate and motion over concatenated trials
    C = diag(C);                                                            % Keep only correlations of each trial with itself
    Corr(i) = mean(C);
end

%% GET MEAN RATE AND MEAN THETA AMPLITUDE PER CELL
mR = cellfun(@(x) nanmean(x,'all'), R);
mAMP = cellfun(@(x) nanmean(x,'all'), AMP);

%% SPLIT MCELLS FROM NON-MCELLS
ismcell = ~isnan(MC(:,1));

%% GROUP MEASURES TO COMPARE
A = {mR, mISI, Corr, mAMP, MDIR/pi, LVEC, ROI};
tit = {'Mean Rate'; 'Median ISI'; 'Speed score'; 'Theta bandpower'; 'Theta phase'; 'Theta modulation'; 'ROI size'};
lA = length(A);

%% COMPARE MEASURES
figure;
for i = 1:lA
    subplot(1,lA,i); hold on
    v = table(A{i}(ismcell)',A{i}(~ismcell)','VariableNames',{'field','no field'});
    violinplot(v);
    p = ranksum(A{i}(ismcell),A{i}(~ismcell));
    if i == 5                                                               % If comparing mean theta phase
        p = circ_wwtest(A{i}(ismcell),A{i}(~ismcell));                      % Use parametric Watson-Williams test .
    end
    plot_significance(p,1,2,A{i}(ismcell),A{i}(~ismcell))
    title(num2str(p));
    ylabel(tit{i});
end
