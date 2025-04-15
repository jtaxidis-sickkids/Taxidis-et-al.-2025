function Figure_field_locomotion(INclass)

%% LOAD POOLED DATA AND STACK
pooldata = ['PooledData_',INclass,'.mat'];
load(pooldata,'R','MOT','MC','time','bins');

% REMOVE MM7-1 5th DAY (NO ODORS) AND MM7-2 3rd DAY (WAS STEP 1)
if strcmp(INclass,'SOM')
    R{1,5} = [];     
    R{2,3} = [];     
    
    MOT{1,5} = [];     
    MOT{2,3} = [];     
    
    MC{1,5} = [];      
    MC{2,3} = []; 
end

R = vertcat(R{:});
MOT = vertcat(MOT{:});
MC = vertcat(MC{:});

ls = length(R);

%% BREAK UP THE DATA CONTAINING MULTIPLE CELLS
count = 1;
Rall = {};
MOTall = {};
MCall = [];

for i = 1:ls
    for r = 1:size(R{i},1)
        Rall{count,1} = squeeze(R{i}(r,:,:))';                          % Keep only that cell's Rates [trials x time]
        MOTall{count,1} = MOT{i};                                           % (Same dimensions)
        MCall(count,:) = MC{i}(r,:);    
        
        count = count + 1;
    end
end
R = Rall;
MOT = MOTall;
MC = MCall;

clear Rall MOTall MCall 

ls = count-1;

%% SMOOTH RATES AND MOTION
for i = 1:ls
    R{i} = smoothdata(R{i},2,'movmean',5);
    MOT{i} = smoothdata(MOT{i},1,'movmean',50);
end

%% KEEP MCELLS ONLY
ismcell = ~isnan(MC(:,1));

R = R(ismcell);
MOT = MOT(ismcell);

ls = sum(ismcell);

%% SPLIT TRIALS INTO LOCOMOTION AND IMMOBILITY DURING ODOR
od1 = (time >=1 & time <= 2);

Rmot = cell(ls,1);
Rimm = cell(ls,1);
Mmot = cell(ls,1);
Mimm = cell(ls,1);
odMOT = cell(ls,1);

for i = 1:ls
    odmot = MOT{i}(od1,:);
    odMOT{i} = nanmean(odmot,1);  
end
figure;histogram(cell2mat(odMOT'),0:0.005:0.3)

for i = 1:ls 
    k = odMOT{i} > 0.02;
    
    Rmot{i} = R{i}(k,:);
    Rimm{i} = R{i}(~k,:);
    
    Mmot{i} = MOT{i}(:,k);
    Mimm{i} = MOT{i}(:,~k);    
end

%% KEEP MEANS PER CELL
Rmot = cellfun(@(x) mean(x,1), Rmot, 'UniformOutput', false);
Rimm = cellfun(@(x) mean(x,1), Rimm, 'UniformOutput', false);

Mmot = cellfun(@(x) mean(x,2), Mmot, 'UniformOutput', false);
Mimm = cellfun(@(x) mean(x,2), Mimm, 'UniformOutput', false);

Rmot = cell2mat(Rmot);
Rimm = cell2mat(Rimm);
Mmot = cell2mat(Mmot')';
Mimm = cell2mat(Mimm')';

%% PLOT 
figure;

for i = 1:ls
    subplot(221); hold on
    plot(time,Mimm(i,:),'Color',[1 0.5 0.5])
    
    subplot(223); hold on
    plot(bins,Rimm(i,:),'Color',[0.8 0.8 0.8])
    
    subplot(222); hold on
    plot(time,Mmot(i,:),'Color',[1 0.5 0.5])
    
    subplot(224); hold on
    plot(bins,Rmot(i,:),'Color',[0.8 0.8 0.8])
end

subplot(221);
fill_plot(time,nanmean(Mimm,1),SEM(Mimm),'r');
xlim([0 3]);

subplot(223);
fill_plot(bins,nanmean(Rimm,1),SEM(Rimm),'k');
xlim([0 3]);

subplot(222);
fill_plot(time,nanmean(Mmot,1),SEM(Mmot),'r');
xlim([0 3]);

subplot(224);
fill_plot(bins,nanmean(Rmot,1),SEM(Rmot),'k');
xlim([0 3]);

