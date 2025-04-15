function Motion_analysis(sessions)

Fs = 1000;
lt = 11*Fs;
t = 0:1/Fs:(lt-1)/Fs;
mV = cell(5,1);

%% LOAD LOCOMOTION
for s = 1:length(sessions)
    s
    dirname = ['../ProcessedData/',sessions{s}];
    load(fullfile(dirname,'CA1data.mat'),'V','opto');

    V = V - mode(V,2);
    V = abs(V);
    
    for i = 0:4
        if i == 0
            k = ~opto(:,1);
        else
            k = opto(:,1) == i & opto(:,2) == 1;
        end
        mV{i+1} = [mV{i+1}; nanmean(V(k,:),1)];      % GET MEAN LOCOMOTION ACROSS TRIALS OF SPECIFIC OPTO
    end
end

%% PLOT
figure;
for i = 1:5
    subplot(5,3,(i-1)*3+1)
    fill_plot(t,nanmean(mV{i},1),SEM(mV{i}),'k');
    axis tight
    xlim([0.5, 9]);
end
