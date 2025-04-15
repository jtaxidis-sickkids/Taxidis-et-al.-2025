session = 'ZD029_0824';

mainfolder = 'Z:\Neuropixels Data Jiannis\Data\PreTraining\Updated_Spike_Sorted';
sessfolder = fullfile(mainfolder,session);
newropixfolder = fullfile(sessfolder,[session,'_g0_imec0']);

dirname = ['../ProcessedData/',session];
load(fullfile(dirname,'ProcessedData.mat'),'pixeldata');
load(fullfile(dirname,'CA1data_new.mat'),'CA1units');

cids = pixeldata.cids;
units = CA1units(:,1);
channels = CA1units(:,2);

Fs = 2500;
time = (-14:25)/Fs * 1000;

for c = 1:length(units)
    k = find(cids == units(c));
    ww = h5read(fullfile(newropixfolder, ['Waveforms_',session,'_Imec0.h5']),sprintf('/%d', k-1));  % Load waveform numbered by unit index in cids (starts from 0)

    channel_vec = channels(c) + 1 + [-8:2:10];
    ww = ww(channel_vec,:,:);
    MW = mean(ww,3);
    
    peak = nan(10,1);
    for ch = 1:10
        peak(ch) = findpeaks(-MW(ch,:),'Npeaks',1,'SortStr','descend');
    end
    [~,maxch] = max(peak);
    
    correct_channel = channel_vec(maxch) - 1;
    
    if maxch ~= 5
        figure('Name',[num2str(units(c)),', ', num2str(channels(c)),', ', num2str(correct_channel)],...
               'units','normalized','outerposition',[0 0 0.2 1]);
        for ch = 1:10
            subplot(10,1,11-ch);hold on;
            for tr = 1:size(ww,3)
                plot(time,ww(ch,:,tr),'Color',[.8 .8 .8]);
            end
            cols =['k';'r'];
            plot(time,MW(ch,:,:),'Color',cols((ch==maxch) + 1),'LineWidth',1.5);
            axis tight
            mygca(ch) = gca;
        end
        yl = cell2mat(get(mygca, 'Ylim'));
        ylnew = [min(yl(:,1)) max(yl(:,2))];
        set(mygca, 'Ylim', ylnew);
    end
end