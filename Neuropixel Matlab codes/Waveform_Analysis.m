function Waveform_Analysis(session)

mainfolder = 'Z:\Neuropixels Data Jiannis\Data\';
sessfolder = fullfile(mainfolder,session);

dirname = ['../ProcessedData/',session];
load(fullfile(dirname,'ProcessedData.mat'),'pixeldata');
load(fullfile(dirname,'CA1data.mat'),'CA1units','chrange','Nc');

if ~exist('Nc','var'), Nc = size(CA1units,1); end

Fs = 2500;
time = (-14:25)/Fs * 1000;


channel_range = [-4:1:5];
lr = length(channel_range);
    
cids = pixeldata.cids;
if ~iscell(cids), cids = {cids}; end
Nprobes = length(cids)

newchannels = [];
W = cell(Nprobes,1);

%%  LOAD WAVEFORMS OF ALL CA1 CELLS AND PLOT PROBLEMATIC ONES
for pr = 1:Nprobes
    newropixfolder = fullfile(sessfolder,[session,'_g0_imec',num2str(pr-1)]);
 
    k = (pr-1)*Nc(1) + (1:Nc(pr));      % ONLY WORKS FOR UP TO 2 PROBES!!
    units = CA1units(k,1);
    channels = CA1units(k,2);
    
    W{pr} = nan(Nc(pr),40,200);                                        % Use it to get number of trials to set size
   
    for c = 1:Nc(pr)
        c
        k = find(cids{pr} == units(c));
        ww = h5read(fullfile(newropixfolder, ['Waveforms_',session,'_Imec',num2str(pr-1),'.h5']),sprintf('/%d', k-1));  % Load waveform numbered by unit index in cids (starts from 0)
       
        nW = size(ww,3);
        channel_vec = channels(c) + 1 + channel_range;
        ww = ww(channel_vec,:,:);
        mww = nanmean(ww,3);
        
        peaks = nan(lr,1);
        for ch = 1:lr
            peaks(ch) = findpeaks(-mww(ch,:),'Npeaks',1,'SortStr','descend');
        end
        [~,maxch] = max(peaks);
        
        proposed_channel = channel_vec(maxch) - 1;
        
        if maxch ~= 5
            fc = figure('Name',[num2str(c),', ', num2str(units(c)),', ', num2str(channels(c)),', ', num2str(proposed_channel)],...
                'units','normalized','outerposition',[0 0 0.2 1]);
            for ch = 1:lr
                subplot(lr/2,2,lr+1-ch);hold on;
                for tr = 1:nW
                    plot(time,ww(ch,:,tr),'Color',[.8 .8 .8]);
                end
                ylabel(['channel ',num2str(ch)])
                cols =['k';'r'];
                h = plot(time,mww(ch,:,:),'Color',cols((ch==maxch) + 1),'LineWidth',1.5);
                if ch == 5, set(h,'Color','b'); end
                axis tight
                mygca(ch) = gca;
            end
            
            yl = cell2mat(get(mygca, 'Ylim'));
            ylnew = [min(yl(:,1)) max(yl(:,2))];
            set(mygca, 'Ylim', ylnew);
            drawnow;
            
            uiwait(fc);
            newchan_indx = input('Enter new channel index: ');
            if isempty(newchan_indx)
                newchan_indx = 5;
            else
                newchan = channel_vec(newchan_indx) - 1
                newchannels = [newchannels; [pr c units(c) channels(c) newchan]];
            end
            
        else
            newchan_indx = 5;
        end
        
        W{pr}(c,:,1:nW) = ww(newchan_indx, :, :);                                      % Keep that channel (+1 cause channels start at 0) x trials x timepoints
    end
end

W = cell2mat(W);
mW = nanmean(W,3);  

%% REMOVE UNITS WITH NEW CHANNELS BEYOND CA1 RANGE
if isempty(newchannels), newchannels = zeros(0,5); end
out = cell(Nprobes,1);
for pr = 1:Nprobes
    k = newchannels(:,1) == pr;
    out{pr} = newchannels(k,5) > chrange(end,pr)+4 | newchannels(k,5) < chrange(1,pr)-4;
end
out = cell2mat(out);
disp(['Number of cells removed = ',num2str(sum(out))]);

%% EXTRACT WAVEFORM FEATURES
Nc = sum(Nc);
K = nan(Nc,2);
Z = nan(Nc,2);
Peak = nan(Nc,1);
Trough = nan(Nc,1);

for c = 1:Nc
    [p1,k1] = findpeaks(-mW(c,:),'Npeaks',1,'SortStr','descend');
    if ~isempty(k1)
        [p2,k2] = findpeaks(mW(c,k1+2:end),'Npeaks',1,'SortStr','descend');
        k2 = k2 + k1+1;
        if ~isempty(k2)
            K(c,:) = [k1 k2];
            
            Peak(c) = p1;                                                 % Waveform peak
            Trough(c) = p2;                                               % Waveform trough
            
            [~,z] = findpeaks(-abs(Peak(c)/2 + mW(c,:)),'SortStr','descend');     % Find times closest to the half-peak points (essentially in their original order)
            Z(c,:) = z(1:2);                                                    % Keep the two closest
        end
    end
end

PTdist = abs(diff(K,[],2))/Fs * 1000;                          % Peak - Trough distance in ms
HW = abs(diff(Z,[],2))/Fs * 1000;                  % Halfwidth in ms

%% PLOT WAVEFORMS
figure;
for c = 1:Nc
    subplot(10,10,c);hold on;
    for tr = 10:10:200
        plot(time,W(c,:,tr),'Color',[.8 .8 .8]);
    end
    plot(time,mW(c,:),'k','LineWidth',1.5);
    title(num2str(CA1units(c,2)));
    drawnow;
    
    if ~isnan(Peak(c))
        plot(time(K(c,1)),-Peak(c),'or','Linewidth',1.5);
        plot(time(K(c,2)),Trough(c),'or','Linewidth',1.5);
        line(time(K(c,1))*[1 1], [0 -Peak(c)],'Color','r','Linewidth',1.5);
        line(time(K(c,2))*[1 1], [0 Trough(c)],'Color','r','Linewidth',1.5);
        line(time([K(c,1) K(c,2)]), -Peak(c)*[1.2 1.2],'Color','k','Linewidth',1.5);
        plot(time([Z(c,1) Z(c,2)]), -Peak(c)/2*[1 1],'Color','g','Linewidth',1.5);
    end
end

%% ORGANIZE WAVEFORM FEATURES
WF = struct('Peak',Peak,'Trough',Trough,'PTdist',PTdist,'HW',HW,...
            'newW',W,'newchannels',newchannels,'out',out);

%% SAVE
dirname = ['../ProcessedData/',session];
save(fullfile(dirname,'CA1data_WF_new.mat'),'WF')


%% JUST PLOT
% dirname = ['../ProcessedData/',session];
% load(fullfile(dirname,'CA1data_WFnew.mat'),'WF')
% Nc = length(WF.Peak);
% figure;
% for c = 1:Nc
%     subplot(10,10,c);hold on;
%     for tr = 10:10:200
%         plot(time,WF.newW(c,:,tr),'Color',[.8 .8 .8]);
%     end
%     plot(time,mean(WF.newW(c,:,:),3),'k','LineWidth',1.5);
% %     title(num2str(WF.newchannels(c)));
%     drawnow;
% 
% %     if ~isnan(WF.Peak(c))
% %         plot(time(K(c,1)),-WF.Peak(c),'or','Linewidth',1.5);
% %         plot(time(K(c,2)),WF.Trough(c),'or','Linewidth',1.5);
% %         line(time(K(c,1))*[1 1], [0 -Peak(c)],'Color','r','Linewidth',1.5);
% %         line(time(K(c,2))*[1 1], [0 Trough(c)],'Color','r','Linewidth',1.5);
% %         line(time([K(c,1) K(c,2)]), -Peak(c)*[1.2 1.2],'Color','k','Linewidth',1.5);
% %         plot(time([Z(c,1) Z(c,2)]), -Peak(c)/2*[1 1],'Color','g','Linewidth',1.5);
% %     end
% end