function CA1ch = Find_CA1(session)

L = matfile(['Z:\Neuropixels Data Jiannis\ProcessedData\',session,'\ProcessedData.mat']);

chV = 2:2:300;                  % Keep every other channel

Fs = 2500;

d60 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',Fs);

d120 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
    'DesignMethod','butter','SampleRate',Fs);

d180 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',179,'HalfPowerFrequency2',181, ...
    'DesignMethod','butter','SampleRate',Fs);


S = L.pixeldata;
LFPs = L.LFPs;

Nprobes = length(S.cgs)

CA1ch = nan(1,Nprobes);

%%
for pr = 1:Nprobes
    %% KEEP GOOD UNITS AND THEIR WAVEFORMS
    goodunits = S.cids{pr}(S.cgs{pr} == 2);                                             % Find cluster IDs with good units
    
    N = length(goodunits);
    channels = zeros(N,1);                                                      % Make matrix for their channel number/depth also stacked by cgs
    
    for i = 1:N                                                                 % For each good unit
        k = find(S.cluster_id{pr} == goodunits(i));                                 % Find where the unit appears in the list of all units
        channels(i) = S.cluster_chan{pr}(k);                                      % Keep the corresponding channel
    end
    
    %% PLOT DISTRIBUTION OF UNITS
    f1 = figure;
    subplot(1,10,1);
    histogram(channels,2:2:300,'orientation','horizontal')
    axis tight
    ylabel('Units per Channel');
    drawnow;
    
    %% LOAD SELECTED LFPS
    lfps = LFPs(chV, 1 : 2000*Fs, pr);
    lfps = double(lfps);
    
    lt = size(lfps,2);
    lfptime = 0 : 1/Fs : (lt-1)/Fs;
    
    %% FILTER LFPS (LOWPASS AND REMOVE 60Hz NOISE)
    for i = 1:length(chV)
        lfp = filtfilt(d60,lfps(i,:));
        lfp = filtfilt(d120,lfp);
        lfps(i,:) = filtfilt(d180,lfp);
    end
    
    %% PLOT SOME LFP SAMPLES
    subplot(1,10,2:8); hold on;
    for i = 1:length(chV)
        plot(lfptime(1:1e5),lfps(i,1:1e5)/max(lfps(i,1:1e5))+i)
    end
    axis tight
    drawnow;
    
    %% DETECT RIPPLES
    [ripples,k,LFPbp] = detect_ripples(lfps,lfptime);
    
    CA1ch(pr) = chV(k)
    
    %% PLOT ALIGNED RIPPLES
    periripdur = 0.1;
    peririp = ceil(periripdur*Fs);
    
    Nch = length(chV);
    rN = size(ripples,1);
    rips = round(ripples*Fs);
    LFPrip = zeros(Nch,2*peririp+1,rN);
    for r = 1:rN                                                                % For each event
        lfptmp = LFPbp(rips(r,1):rips(r,2));                                    % Keep the corresponding bandpassed signal
        [~,rm] = max(lfptmp);                                                   % Find the position of its maximum value
        rmax = rips(r,1) - 1 + rm;                                           % Store the overall array position of the peak
        
        LFPrip(:,:,r) = lfps(:,rmax + (-peririp : peririp));
    end
    
    figure(f1);
    subplot(1,10,9:10); hold on;
    perit = (-peririp:peririp)/Fs;
    step = 10;
    for ch = 1:Nch
        lfp = squeeze(LFPrip(ch,:,:))';
        plot(perit,mean(lfp) - (mean(lfp(:,1))) + step*ch);
        if chV(ch) == CA1ch(pr)
            plot(perit,mean(lfp) - (mean(lfp(:,1))) +step*ch,'k','LineWidth',2);% SEM(lfp),'k')
        end
    end
    axis tight;
    ylabel('Pyramidal Layer');
    set(gca,'Ytick',step:step:step*length(chV),'YTickLabel',string(chV));
    drawnow;
end