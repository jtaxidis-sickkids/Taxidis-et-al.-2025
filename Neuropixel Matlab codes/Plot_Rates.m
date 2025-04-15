function Plot_Rates(session)

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

%% GET FIRST ODOR
firstodor = floor(progress(:,1)/10);                                        % Keep each trial's first odor
Ntr = length(firstodor);

%% COMPUTE MEAN RATES AND LOCOMOTION
mR = mean(R,3);                                                  % Get mean over trials
sR = std(R,[],3) / sqrt(Ntr);                              % and SEM

mV = mean(V,1);
sV = SEM(V);

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
            sp = CA1spikes{cc,tr};
            plot(sp,ones(1,length(sp))+tr-1,'.','Markersize',3,'Color',C(firstodor(tr),:))
        end
        xlim([0 11]);
        ylim([0 Ntr]);
        ylabel([num2str(PY_IN(cc)),',',num2str(CA1units(cc,2)),',',num2str(CA1units(cc,3))]);

        % Plot average rates
        subplot(10,2,2*c); hold on;
        fill_plot(bins(1:end-1),mR(cc,:),sR(cc,:),[0.8 0.8 0.8]);

        xlim([0 11]);
    end
    drawnow;
end

