function Organize_SpikeSignal(days)

Fs = 30.9;                                                                  % Scope samplign rate

for dd = 1:length(days)                                                     % For each day
    load([days{dd},'/Behavior.mat'],'Params','Motion','trial_on');
    load([days{dd},'/Rates.mat'],'Spikes','PY');
    
    ls = length(Spikes);                                                     % Number of sessions
    Ns = size(Spikes{1},1);                                                  % Number of segments (incluing PY and IN)
    Ntr = cellfun(@length, trial_on);                                       % Number of trials (BASED ON EDR)
    
    delay = zeros(1,ls);
    for s = 1:ls                                                            % For each session
        if ~isempty(Params{s})                                              % If that session was not missing a Mat file
            delay(s) = Params{s}.delay;                                     % Get the delay in that session
        else
            delay(s) = nan;
        end
    end
    delays = unique(delay);                                                 % Keep the different delay values used
    delays(isnan(delays)) = [];                                             % Remove nans from trials with missing MAT files
    ld = length(delays);                                                    % Number of different delays
    
    for d = 1:ld                                                            % For each delay
        %% GET SESSIONS OF SAME DELAY
        disp(['Delay = ',num2str(delays(d))]);
        
        ses_del = find(delay == delays(d));                                 % Find which sessions had that delay
        lsd = length(ses_del);                                              % and count them
        
        mot = Motion(ses_del);                                              % and motion
        spikes = Spikes(ses_del);                                           % deconvolved (spiking) signal
        tr_on = trial_on(ses_del);                                          % and trial onset points (for DFF and Rates)
        
        % THIS DAY, WINEDR WAS STOPPED BEFORE VIDEO. ADD ZEROS AT THE END OF FINAL-SESSION MOTION
        if strcmp(days{dd} , 'I:\Jiannis - Imaging Analysis\JT87\Processed_files\2017_12_23')
            mot{6} = [mot{6}; zeros(size(spikes{6},2) - length(mot{6}),1)];
        end
        
        %% FIND MINIMUM TRIAL DURATION (EXCLUDING VIDEOS THAT WERE CUT SHORT)
        tr_on_temp = cell(size(tr_on));
        for s = 1:lsd                                                       % For each session of that delay
            tr_on_temp{s} = [tr_on{s}(:,1); length(mot{s})];                % Concatenate DFF trial onset points, and session lengths
        end
        dur = cellfun(@diff, tr_on_temp,'UniformOutput',false);             % Compute duration of each trial
        dur = cell2mat(dur);
        
        [sorted_dur,order] = sort(dur);                                     % Sort the DFF datapoint durations of all trials
        trial_datapoints = Params{ses_del(1)}.trial_dur * Fs;               % Get the datapoints of the actual trial of that delay
        k = find(sorted_dur > trial_datapoints,1,'first');                  % Find the minimum EDR trial duration that is longer than the behavioral trial duration
        dur = dur(order);                                                   % (which excludes videos that were accidentally stopped before the final trial was over)
        mindur = dur(k);                                                    % Keep minimum durations for DFF
        
        %% ADD ZEROS IF SOME FINAL TRIAL WAS CUT SHORT
        for s = 1:lsd                                                       % For each session of that delay
            ntr = Ntr(ses_del(s));                                          % Keep the number of trials
            endpoint = tr_on{s}(ntr,1) + mindur -1;                       % Find what the final datapoint shoudl be according to the minimum duration
            if size(spikes{s},2) < endpoint                                    % If the DFF was cut short 
                spikes{s} = [spikes{s}, zeros(Ns, endpoint-size(spikes{s},2))];% add zeros at the end                
            end
        end
    
        %% ORGANIZE BY TRIALS
        S = cell(1,1,lsd);
        for s = 1:lsd                                                       % For each session of that delay
            S{s} = zeros(Ns,mindur,Ntr(ses_del(s)));                        % Allocate memory for spike signal (Ns x time x Ntr)
            
            for tr = 1:Ntr(ses_del(s))                                      % For each trial in that session
                S{s}(:,:,tr) = spikes{s}(:,tr_on{s}(tr,1) + (0:mindur-1));  % Keep the spike singal of that trial
            end
        end
        S = cell2mat(S);                                                    % Concatenate all sessions together
        
        %% REARRANGE BY CELL TYPE
        S = [S(PY == 1,:,:); S(PY == 0,:,:)];                               % Keep ALL PY on top with all IN underneath      
        S = real(S);

        %% SAVE DATA FOR EACH DELAY
        [animal,main_folder,day] = get_animal_name(days{dd});
        dirname = fullfile(main_folder,[animal,'_Analysis'],day);           % Create directory
        if ~exist(dirname,'dir')                                            % With 'mouse name'_Analysis
            mkdir(dirname);
        end
        save(fullfile(dirname,['Trials_delay_',num2str(delays(d)),'_SpikeSignal.mat']),'S');        % Save data from all sessions with this delay
    end
end

