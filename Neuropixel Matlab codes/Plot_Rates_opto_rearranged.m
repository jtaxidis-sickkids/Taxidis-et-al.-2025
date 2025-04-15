function Plot_Rates_opto_rearranged(session)

dirname = ['../ProcessedData/',session];
load(fullfile(dirname,'CA1data.mat'),'CA1spikes','R','opto','progress','V','bins','CA1units');
load(fullfile(dirname,'PY_IN.mat'),'PY_IN'); 
load(fullfile(dirname,'Mcells.mat'),'Mcells');

C = cell_colors;
C = [C(1:2,:); .4 .4 .4; C(3,:)];

[Nc,Ntr] = size(CA1spikes);
Nc
sum(PY_IN == 1)
sum(PY_IN == 2)
sum(isnan(PY_IN))

%% ARRANGE CELLS BY CELLTYPE
R = cat(1, R(PY_IN == 1,:,:), R(PY_IN == 2,:,:), R(isnan(PY_IN),:,:)); 
CA1spikes = cat(1, CA1spikes(PY_IN == 1,:), CA1spikes(PY_IN == 2,:), CA1spikes(isnan(PY_IN),:)); 
CA1units = cat(1, CA1units(PY_IN == 1,:), CA1units(PY_IN == 2,:), CA1units(isnan(PY_IN),:)); 

for tt = 1:2
    if isempty(Mcells{2,tt}), Mcells{2,tt} = zeros(0,3); end
    Mcells{2,tt}(:,1) = Mcells{2,tt}(:,1) + sum(PY_IN==1);
end

PY_IN = sort(PY_IN); % nans go to bottom

%% SMOOTH RATES
R = smoothdata(R,2,'gaussian',12); 

%% ARRANGE TRIALS BY NO LIGHT AND THEN BY OPTO PROTOCOL (GROUPS TOGETHER ALL REPETITIONS OF SAME PROTOCOL)
opto(~opto(:,2), 1) = 0;                % Set all trials without light as opto = 0
% opto(opto(:,1) == 4, 1) = 0.5;          % Set all trials with full-odor light as opto = 0.5
[~,order] = sort(opto(:,1));            % Then sort by opto

opto = opto(order,:);
CA1spikes = CA1spikes(:,order);
R = R(:,:,order);
V = V(order,:);
progress = progress(order,:);

prot = opto(:,1);
[protunique,init_trials] = unique(prot); 
init_trials = [init_trials; Ntr+1];
lprot = length(protunique);

%% GET FIRST ODOR AND ON/OFF CHANGE TRIALS
firstodor = floor(progress(:,1)/10);

%% COMPUTE MEAN RATES AND LOCOMOTION
mR = cell(lprot,1);
sR = cell(lprot,1);
mV = cell(lprot,1);
sV = cell(lprot,1);

for i = 1:lprot
    k = (prot == protunique(i));
    Rtemp = R(:,:,k);
    Vtemp = V(k,:);
    
    mR{i} = mean(Rtemp,3);
    sR{i} = std(Rtemp,[],3)/sqrt(sum(k));
    
    mV{i} = mean(Vtemp,1);
    sV{i} = SEM(Vtemp);
    
    if isempty(Rtemp), mR{i} = [];  sR{i} = [];  mV{i} = [];  sV{i} = []; end
end

%% PLOT SPIKES AND MEAN RATES PER CELL AND LOCOMOTION
for f = 1:ceil(Nc/10)
    figure;
    cmax  = 10;
    if f == ceil(Nc/10), cmax = Nc - (f-1)*10; end
    
    for c = 1:cmax
        % Plot spikes
        subplot(10,2*lprot-1,(2*lprot-1)*(c-1) + (1:lprot)); hold on;
        line([1 1],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)                                      
        line([2 2],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)
        line([7 7],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)
        line([8 8],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)

        for i = 2:lprot
            line([0 11],init_trials(i)*[1 1],'Color','k','LineWidth',1);
        end
 
        cc  = (f-1)*10 + c;
        for tr = 1:Ntr
            col = C(firstodor(tr),:);
            sp = CA1spikes{cc,tr};
            plot(sp,ones(1,length(sp))+tr-1,'.','Markersize',6,'Color',col)
        end
        xlim([0 11]);
        ylim([0 Ntr+1]);
        ylabel([num2str(PY_IN(cc)),',',num2str(CA1units(cc,2)),',',num2str(CA1units(cc,3))]);

        % Plot opto protocol 
        for i = 2:lprot
            X = (protunique(i) == 1)*([1.03 1.03 1.05 1.05])...
                + (protunique(i) == 2)*([1.05 1.05 1.15 1.15])...
                + (protunique(i) == 3)*([1.2 1.2 1.4 1.4])...
                + (protunique(i) == 4)*([1 1 2 2]);       
            Y = [init_trials(i) init_trials(i+1) init_trials(i+1) init_trials(i)];
            patch(X,Y,'b','EdgeColor','none', 'FaceAlpha',0.2);
        end
        
         % Plot field
         if size(Mcells,2) == 2         % IF PLOTTING 2Pversion MCELLS
             if ~isnan(PY_IN(cc))
                 for tt = 1:2
                     k = find(Mcells{PY_IN(cc),tt}(:,1) == cc);
                     if ~isempty(k)
                         line(Mcells{PY_IN(cc),tt}(k,2)* [1 1], [0 Ntr],'Color',C(tt,:),'LineWidth',1.5);
                     end
                 end
             end
         end
         
        % Plot rates with and without opto for each protorol
        subplot(10,(2*lprot-1),(2*lprot-1)*(c-1) + lprot + 1); hold on;
        fill_plot(bins(1:end-1),mR{1}(cc,:),sR{1}(cc,:),C(3,:));
        fill_plot(bins(1:end-1),mR{4}(cc,:),sR{4}(cc,:),C(4,:));
        xlim([0 11]);

        for i = 2:lprot-1
            subplot(10,2*lprot-1,(2*lprot-1)*(c-1) + lprot + i); hold on;
            fill_plot(bins(1:end-1),mR{i}(cc,:),sR{i}(cc,:),C(4,:));
        end
        xlim([0 11]);
    end
    drawnow; 
end
 
return

%% PLOT MEAN RATES
figure;
subplot(lprot-1,1,1); hold on;
for c = 1:Nc
    plot(bins(1:end-1),mR{1}(c,:),'Color',C(3,:));
    plot(bins(1:end-1),mR{4}(c,:),'Color',C(4,:));
end
fill_plot(bins(1:end-1),mean(mR{1},1),SEM(mR{1}),'k');
fill_plot(bins(1:end-1),mean(mR{4},1),SEM(mR{4}),'b');
xlim([0 11]);

for i = 2:lprot-1
    for c = 1:Nc
        subplot(lprot-1,1,i); hold on;
        plot(bins(1:end-1),mR{i}(c,:),'Color',C(4,:));
    end
    fill_plot(bins(1:end-1),mean(mR{i},1),SEM(mR{i}),'b');
    xlim([0 11]);
end

%% PLOT MEAN LOCOMOTION
figure;
subplot(lprot-1,1,1); hold on;
fill_plot(0:0.001:11-0.001,mV{1},sV{1},'k');
fill_plot(0:0.001:11-0.001,mV{4},sV{4},'b');
xlim([0 11]);

for i = 2:lprot-1
    subplot(lprot-1,1,i); hold on;
    fill_plot(0:0.001:11-0.001,mV{i},sV{i},'b');
    xlim([0 11]);
end




