function VOLPY_processing_BASIC(GEVI,ID,day,session,Fs)

%% LOAD VOLPY DATA
[asapfile,asapfiles,video] = get_GEVIfile(GEVI,ID,day,session);
disp(asapfile);

[path,videoname] = fileparts(asapfile);
VOLPYfile = fullfile(path,'Volpy_data',[videoname,'_volpy_data.mat']);
load(VOLPYfile,'DFF','ROIs','spikes','weights','mean_im');

[~,Nfr] = size(DFF);

[Nr,Nfr] = size(DFF);

%% SET TRIAL DURATION
trial_dur = 30; 
if contains(videoname,'long')
    trial_dur = 16;
end

ls = trial_dur*Fs;                         % Length of section (trial)
ASAPtime = (0:ls-1)/Fs;

Ntr = ceil(Nfr/ls);                             % Number of trials
disp(['Number of trials = ',num2str(Ntr)]);

%% RENAME AND RESHAPE DATA
if iscell(mean_im)
    mean_im = mean_im{1};
    mean_im = reshape(mean_im,[1,size(mean_im)]);
end
Template = squeeze(mean_im(1,:,:));

CC = cell(Nr,1);
for r = 1:Nr
    structBoundaries = bwboundaries(squeeze(ROIs(r,:,:)));                  % Get coordinates of the boundary of the freehand drawn region.
    k = structBoundaries{1};                                                % Get n by 2 array of x,y coordinates.
    CC{r} = flipud(k');
end    

SP = cell(Nr,Ntr);
if ~iscell(spikes)
    sptemp = spikes;
    spikes = {};
    for r = 1:Nr
        spikes{r} = sptemp(r,:);
    end
end

totspikes = 0;
for r = 1:Nr
    for tr = 1:Ntr
        k = spikes{r} > (tr-1)*ls & spikes{r} <= tr*ls;        
        sp_tr = spikes{r}(k) - (tr-1)*ls;
        SP{r,tr} = ASAPtime(sp_tr);
        totspikes = totspikes + length(SP{r,tr});
    end
    SP{r,1}(SP{r,1} < 0.2) = [];                                            % REMOVE SPIKES BEFORE THE FIRST 200ms FROM FIRST TRIAL
end
totspikes

if Nfr < Ntr*ls
    DFF = cat(2,DFF,nan(Nr, Ntr*ls - Nfr));
end
DFF = reshape(DFF,[Nr,ls,Ntr]);
DFF = DFF * 100;                                                            % Turn to percentage

%% COMPUTE FIRING RATES
binl = 0.1;                                                                 % Firing rate time bin 100ms
bins = 0 : binl : trial_dur;

R = nan(Nr,length(bins)-1,Ntr);
for r = 1:Nr
    for tr = 1:Ntr
        R(r,:,tr) = histcounts(SP{r,tr},bins) / binl;                       % Compute firing rate
%         R(r,:,tr) = smooth(R(r,:,tr),7);                                  % And smooth it
    end
end

bins = bins + binl/2;
bins(end) = [];


%% ORGANIZE FINAL STRUCTURES
V.Template = Template;
V.weights = weights;   
V.ROI = ROIs;
V.CC = CC;
V.DFF = DFF;
V.time = ASAPtime;
V.SP = SP;
V.R = R;
V.bins = bins;

%% STORE DATA 
outfile = fullfile(path,[videoname,'_data.mat']);
save(outfile,'V'); 

%% PLOT PROCESSING
% Plot_VOLPY_processing(ID,day,session);
Template = V.Template;
weights = V.weights;   
ROI = V.ROI;
CC = V.CC;
DFF = V.DFF;
ASAPtime = V.time;
SP = V.SP;

%% CREATE FIGURE
fig = figure('Name',[ID,', ',day,', session ',num2str(session),' (',videoname,')']);
screensize = get(0,'Screensize');
set(gcf, 'PaperUnits', 'points', 'Units', 'points');
set(gcf, 'Position', [10 0 screensize(3)*0.9 screensize(4)*0.7]);

%% PLOT MOTION CORRECTION
colormap gray;

[Nr,~,Ntr] = size(DFF);
lr = max(6+2*Nr,Ntr);

subplot(lr,10,[1:4,11:14,21:24]);
imshow(Template, [min(Template(:)), max(Template(:))]);
axis equal; title('Motion-corrected');

%% PLOT ROI
for r = 1:Nr
    subplot(lr,10,[2*r*10+(11:14) , 2*r*10+(21:24)]);
    imshow(squeeze(weights(r,:,:)),[min(weights(:)), max(weights(:))]);
    axis equal; title('Final ROI');
    hold on;
    
    plot(CC{r}(1,:),CC{r}(2,:),'-w');
    C = regionprops(squeeze(ROI(r,:,:)), 'Centroid');
    text(C.Centroid(1),C.Centroid(2),num2str(r),'Color','r','Fontsize',16);
end

%% PLOT TRACES AND SPIKES PER TRIAL
C = [255 195 55;    % yellow
    36 146 102;     % green
    9 117 220]/255;  % red

Cp = parula(ceil(1.5*Nr));
    
ymin = -1.1; 
ymax = (Nr-1)*2.5 + 1.5;
yshade = [ymin ymax ymax ymin];

for tr = 1:Ntr
    % PLOT TRIAL ------------------
    subplot(lr,10, (tr-1)*10 + (5:9)); hold on;
    for r = 1:Nr
        plot(ASAPtime, DFF(r,:,tr)  + (r-1)*2.5,'Color',Cp(r,:));
        plot(SP{r,tr},zeros(size(SP{r,tr})) + (r-1)*2.5 + 1.2,'.k');
    end
    % ----------------------
    
    % PLOT PERISPIKE ----------
    subplot(lr,10, tr*10); hold on;
    window = 20;
    peritime = -window/Fs : 1/Fs : window/Fs;
    for r = 1:Nr
        perispike = nan(length(SP{r,tr}),2*window+1);
        for sp = 1:length(SP{r,tr})
            spiketime = SP{r,tr}(sp);
            spikeindex = find(ASAPtime == spiketime);
            spikewindow = spikeindex + (-window:window);
            if spikewindow(1) >= 1 && spikewindow(end) <= length(ASAPtime)
                perispike(sp,:) = DFF(r,spikewindow,tr);
%                 perispike(sp,:) = perispike(sp,:) - mean(perispike(sp,1:5));   % REMOVE THE FIRST 5ms TO SET AS BASELINE
            end
            plot(peritime,perispike(sp,:) + 2*(r-1),'Color',[0.5 0.5 0.5],'Linewidth',0.5);
        end
        if ~isempty(perispike)
            plot(peritime,nanmean(perispike,1) + 2*(r-1),'Color',Cp(r,:),'Linewidth',3);
            axis tight;
        end
    end
end
drawnow;

