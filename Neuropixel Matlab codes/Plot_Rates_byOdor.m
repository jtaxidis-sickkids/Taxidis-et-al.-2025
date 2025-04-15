function Plot_Rates_byOdor(session)

dirname = ['../ProcessedData/',session];
load(fullfile(dirname,'CA1data.mat'),'CA1spikes','R','progress','V','bins','CA1units'); 
load(fullfile(dirname,'PY_IN.mat'),'PY_IN'); 
load(fullfile(dirname,'Mcells.mat'),'Mcells');

C = cell_colors;
C = C(1:2,:);

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

%% ARRANGE BY FIRST ODOR
firstodor = floor(progress(:,1)/10);                                        % Keep each trial's first odor

[fod,order] = sort(firstodor);                                              % Sort the first odor
lastA = sum(fod == 1);                                                      % Count odorA trials

CA1spikes = CA1spikes(:,order);                                             % Rearrange spike trials according to 1st odor
R = R(:,:,order);                                                           % and rates
V = V(order,:);                                                             % and locomotion

%% COMPUTE MEAN RATES AND LOCOMOTION
mR = cell(2,1);
sR = cell(2,1);
mV = cell(2,1);
sV = cell(2,1);

for i = 1:2                                                                 % For each 1st odor
    k = (1:lastA)*(i==1) + ((lastA+1):Ntr)*(i==2);                          % keep corresponding trial indexes
    
    Rtemp = R(:,:,k);                                                       % Keep corresponding rates
    mR{i} = mean(Rtemp,3);                                                  % Get mean over trials
    sR{i} = std(Rtemp,[],3) / sqrt(length(k));                              % and SEM
    
    Vtemp = V(k,:);                                                         % Same for locomotion
    mV{i} = mean(Vtemp,1);
    sV{i} = SEM(Vtemp);
end

%% PLOT SPIKES AND MEAN RATES PER CELL AND LOCOMOTION
for f = 1:ceil(Nc/10)
    figure;
    cmax  = 10;
    if f == ceil(Nc/10), cmax = Nc - (f-1)*10; end
    
    for c = 1:cmax
        % Plot spikes
        subplot(10,2,2*c-1); hold on;
        line([1 1],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)
        line([2 2],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)
        line([7 7],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)
        line([8 8],[0 Ntr+1],'Color',[.8 .8 .8],'LineWidth',1)
        
        cc  = (f-1)*10 + c;                 % Actual cell to plot
        for tr = 1:Ntr
            col = C((tr > lastA) + 1,:);
            sp = CA1spikes{cc,tr};
            plot(sp,ones(1,length(sp))+tr-1,'.','Markersize',3,'Color',col)
        end
        xlim([0 11]);
        ylim([0 Ntr]);
        ylabel([num2str(PY_IN(cc)),',',num2str(CA1units(cc,2)),',',num2str(CA1units(cc,3))]);
        
        % Plot field
        if ~isnan(PY_IN(cc))
            for tt = 1:2
                k = find(Mcells{PY_IN(cc),tt}(:,1) == cc);
                if ~isempty(k)
                    line(Mcells{PY_IN(cc),tt}(k,2)* [1 1], [0 lastA]*(tt==1) + [lastA+1 Ntr]*(tt==2),'Color','k','LineWidth',1.5);
                end
            end
        end
        
        % Plot average rates
        subplot(10,2,2*c); hold on;
        for i = 1:2
            fill_plot(bins(1:end-1),mR{i}(cc,:),sR{i}(cc,:),C(i,:));
        end
        xlim([0 11]);
    end
    drawnow;
end

%% PLOT POOLED MEAN RATES AND LOCOMOTION
figure;
for i = 1:2
    for c = 1:Nc
        subplot(3,1,i); hold on;
        plot(bins(1:end-1),mR{i}(c,:),'Color',[.8 .8 .8]);
    end
    
    subplot(3,1,i);
    fill_plot(bins(1:end-1),mean(mR{i},1),SEM(mR{i}),C(i,:));
    xlim([0 11]);
    
    subplot(313); hold on
    fill_plot(0:0.001:11-0.001,mV{i},sV{i},C(i,:));
    xlim([0 11]);
end

