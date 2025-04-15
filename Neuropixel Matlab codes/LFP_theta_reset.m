function LFP_theta_reset(sessions)

addpath('CircStat2012a');

%% MAKE FILTERS
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

[filtb,filta] = butter(3,[4 10]/(Fs/2),'bandpass');                        % construct the bandpass filter

lt = 11*Fs;
t = 0:1/Fs:(lt-1)/Fs;

mL = cell(5,1);
for i = 1:5
    mL{i} = nan(length(sessions),lt);
end
mfL = mL;
vPH = mL;

counter = 1;

%% LOAD LFPS AND FILTER
for s = 1:length(sessions)
    s
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'L','opto','chrange','CA1ch');
    
    [~,~,Ntr,lprobes] = size(L);
    
    for j = 1:lprobes
        Ltemp = L(chrange(:,j) == CA1ch(j),:,:,j);  % Keep the CA1 pyr layer LFP
        Ltemp = squeeze(Ltemp);        % Turn to time x trials
        Ltemp = double(Ltemp);
        
        % FILTER OUT 60Hz NOISE
        Ltemp = Ltemp(:)';  % Concatenate trials and turn to single row
        Ltemp = filtfilt(d60,Ltemp);
        Ltemp = filtfilt(d120,Ltemp);
        Ltemp = filtfilt(d180,Ltemp);
        
        % THETA BANDPASS
        fL = filtfilt(filtb,filta,Ltemp);
        
        H = hilbert(fL);                                                  % Get hilbert transform
        phase = angle(H);                                            % Get phases
        
        % RESHAPE INTO TRIALS
        Ltemp = reshape(Ltemp,lt,Ntr);              % Reshape back to time x trials
        fL = reshape(fL,lt,Ntr);              % Reshape back to time x trials
        phase = reshape(phase,lt,Ntr);              % Reshape back to time x trials
        phase = real(phase);
        
        for i = 0:4
            if i == 0
                k = ~opto(:,1);
            else
                k = opto(:,1) == i & opto(:,2) == 1;
            end
            
            mL{i+1}(counter,:) = nanmean(Ltemp(:,k),2);      % GET MEAN LFP ACROSS TRIALS OF SPECIFIC OPTO
            mfL{i+1}(counter,:) = nanmean(fL(:,k),2);      % GET MEAN LFP ACROSS TRIALS OF SPECIFIC OPTO
            vPH{i+1}(counter,:) = circ_var(phase(:,k),[],[],2);      % GET THETA PHASE VARIANCE ACROSS TRIALS OF SPECIFIC OPTO
        end
        
        counter = counter + 1;
    end
end

%% PLOT
figure;
for i = 1:5
    subplot(5,3,(i-1)*3+1)
    fill_plot(t,nanmean(mL{i},1),SEM(mL{i}),'k');
    axis tight
    xlim([0.5, 9]);
    
    subplot(5,3,(i-1)*3+2)
    fill_plot(t,nanmean(mfL{i},1),SEM(mfL{i}),'k');
    axis tight
    xlim([0.5, 9]);
    
    subplot(5,3,i*3)
    fill_plot(t,nanmean(vPH{i},1),SEM(vPH{i}),'k');
    axis tight
    xlim([0.5, 9]);
end

%%
% fvec = 0:0.1:100;
% % COMPUTE PSD
% [~,~,tf,psd] = spectrogram(L, 512, 256, fvec, Fs,'yaxis');
% % Interpolate
% tq = tf(1):0.1:tf(end);
% fq = fvec(1):0.025:fvec(end);
% [T,Fint] = meshgrid(tf,fvec);
% [Tq,Fq] = meshgrid(tq,fq);
% psdq = interp2(T,Fint,psd,Tq,Fq);
% % Smooth
% psdq = imgaussfilt(psdq,'FilterSize',5);
% PSD = psdq;
% imagesc(tq,fq,PSD,[min(PSD(:)),max(PSD(:))])
% -------------

% COMPUTE SPECTRUM
% [Slfp,Fvec] = pspectrum(L,Fs,'FrequencyLimits',[0 200],'FrequencyResolution',0.5);      % Compute power spectrum during motion
% Slfp(~isnan(Slfp)) = pow2db(Slfp(~isnan(Slfp)));                                    % Switch to decibel (skip Nan entries where it crashes)

