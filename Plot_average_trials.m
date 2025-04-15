function Plot_average_trials(ID,day,session)

%% LOAD DATA
[asapfile,asapfiles] = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_data.mat']);
load(Datafile);

pars = B.pars;
bins = V.bins;
ASAPtime = V.time;

% LOAD THETA-MOTION (ALREADY POOLED ANALYSED)
datafile = fullfile(path,[videoname,'_Freq_Motion_data.mat']);
load(datafile,'F');

%% POOL TRIALS AND KEEP UP TO MAX TRIAL
[DFF,SP,R,MOT,trials,outcomes] = load_and_pool(ID,day,session);                               % Pool trials of the same ROI

[Nr,Nt,Ntr] = size(DFF);
lf = size(F.Ampl,3);

mtr = max(trials);

%% INITIALIZE
EDR_trial_dur = 11 * 1e3;                                                   % EDR-trial duration (Fs = 1e3)
EDRtime = 0 : 1e-3 : (EDR_trial_dur-1)*1e-3;

C = [9 117 220;  % blue
    255 195 55;    % yellow
    36 146 102;     % green
    200 150 140]/255; % red

corr_err = zeros(1,Ntr);
corr_err(outcomes == 1 | outcomes == 4) = 1;
corr_err(outcomes == 2 | outcomes == 3) = 2;
outcome_tit = {'Correct'; 'Error'};

tmpts = [pars.stim1_on, pars.delay_on, pars.stim2_on, pars.stim2_off, ...
    pars.lick_on, pars.vacuum_on,  pars.trial_dur];

bins(end) = [];

%% RESHAPE BANDPASSED VARIABLES
Ampl = reshape(F.Ampl,[Nr,Nt,Ntr,lf]);            % Split trials into the third dimension
Phase = reshape(F.Phase,[Nr,Nt,Ntr,lf]);            % Split trials into the third dimension

%% DOWNSAMPLE TO 500Hz
MOT = MOT(1:2:end,:);
DFF = DFF(:,1:2:end,:);
Ampl = Ampl(:,1:2:end,:,:);
Phase = Phase(:,1:2:end,:,:);

ASAPtime = ASAPtime(1:2:end);
EDRtime = EDRtime(1:2:end);

%% NORMALIZE
MOT = (MOT-min(MOT(:))) / max(MOT(:));      
a1 = find(ASAPtime >= 1,1,'first');      
a2 = find(ASAPtime <= 10,1,'last');
b1 = find(bins >= 1,1,'first');      
b2 = find(bins <= 10,1,'last');

for r = 1:Nr
    DFF(r,:,:) = (DFF(r,:,:)-min(abs(DFF(r,a1:a2,:)),[],'all')) / max(abs(DFF(r,a1:a2,:)),[],'all');
    R(r,:,:) = R(r,:,:) / max(R(r,b1:b2,:),[],'all');           % R already has 0 as minimum
    for f = 1:lf
        Ampl(r,:,:,f) = (Ampl(r,:,:,f) - min(Ampl(r,a1:a2,:,f),[],'all')) / max(Ampl(r,a1:a2,:,f),[],'all');
        Phase(r,:,:,f) = (Phase(r,:,:,f) - min(abs(Phase(r,:,:,f)),[],'all')) / max(abs(Phase(r,:,:,f)),[],'all') / 5;
    end
end

%% PLOT ALL TRIALS
for r = 1:Nr
    figure('Name',[ID,', ',day,', session ',num2str(session),' ROI ',num2str(r)]);
    
    for i = 0:mtr
        subplot(1,mtr+1,i+1); hold on;
        if i == 0                   % PLOT ALL TRIALS
            title('All Trials');
            tt = 1:length(trials);
        else                        % PLOT BY TRIAL TYPE
            tt = find(trials == i);
            title(['ODOR-',num2str(i)]);
        end
        Ntt = length(tt);      
        
        ymax = lf+4 + (Ntt+1)*0.05 ;
        yshade = [0 ymax ymax 0];
        patch([tmpts(1)*[1 1] tmpts(2)*[1 1]], yshade, C(i+1,:),'EdgeColor','none', 'FaceAlpha',0.2);
        patch([tmpts(3)*[1 1] tmpts(4)*[1 1]], yshade, C(1,:),'EdgeColor','none', 'FaceAlpha',0.2);
        
        % PLOT PER TRIAL
        for tr = 1:Ntt
            % PLOT DFFs, THETA AND RATES
            plot(EDRtime',MOT(:,tt(tr)),'Color',C(4,:),'LineWidth',0.1);
            plot(ASAPtime,DFF(r,:,tt(tr)) + 2, 'Color', C(1,:),'LineWidth',0.1);
            plot(bins,R(r,:,tt(tr)) + lf+3, 'Color',[0.8 0.8 0.8],'LineWidth',0.1);
            for f = 1:lf
%                 plot(ASAPtime,Phase(r,:,tt(tr),f) + 2*f+1.5, 'Color', C(1,:),'LineWidth',0.1);
                plot(ASAPtime,Ampl(r,:,tt(tr),f) + f+2, 'Color', [0.8 0.8 0.8],'LineWidth',0.1);
            end
            % PLOT SPIKES
            sp = SP{r,tt(tr)};
            if ~isempty(sp)
                plot(sp,zeros(size(sp)) + lf+4 + tr*0.05, '.','Color',[.5 .5 .5],'Markersize',4);
            end
        end
        
        % PLOT MEANS
        fill_plot(EDRtime',mean(MOT(:,tt),2),SEM(MOT(:,tt)),'r');
        fill_plot(ASAPtime,mean(DFF(r,:,tt),3) + 2, SEM(squeeze(DFF(r,:,tt))'), 'b');
        fill_plot(bins,mean(R(r,:,tt),3) + lf+3, SEM(squeeze(R(r,:,tt))'), 'k');
        for f = 1:lf
%             fill_plot(ASAPtime,mean(Phase(r,:,tt,f),3) + 2*f+1.5, SEM(squeeze(Phase(r,:,tt,f))'), 'k');
            fill_plot(ASAPtime,mean(Ampl(r,:,tt,f),3) + f+2, SEM(squeeze(Ampl(r,:,tt,f))'), 'k');
        end
        
        xlim([EDRtime(1) EDRtime(end)]); 
        ylim([0 ymax]);
    end
end

%% PLOT MEAN TRIALS BY OUTCOME
% for i = 1:2
%     subplot(2,2,2*i); hold on;
%     patch([tmpts(1)*[1 1] tmpts(2)*[1 1]], yshade, C(3,:),'EdgeColor','none', 'FaceAlpha',0.2);
%     patch([tmpts(3)*[1 1] tmpts(4)*[1 1]], yshade, C(3,:),'EdgeColor','none', 'FaceAlpha',0.2);
%
%     k = find(corr_err == i);
%     le = sqrt(length(k));
%
%     fill_plot(EDRtime',mean(MOT(:,k),2),std(MOT(:,k),[],2)/le ,'r');
%     for r = 1:Nr
%         fill_plot(ASAPtime,mean(DFF(r,:,k),3)+ r-1/2, std(DFF(r,:,k),[],3)/le , Cp(r,:));
% %         sp = cell2mat(SP(r,k));
% %         if ~isempty(sp)
% %             plot(sp,zeros(size(sp))+ r,'.k');
% %         end
%         fill_plot(bins,mean(R(r,:,k),3)+r+0.1, std(R(r,:,k),[],3)/le, 'k');
%     end
%     xlim([EDRtime(1) EDRtime(end)]); ylim([ymin ymax]);
%     title(outcome_tit{i});
% end
