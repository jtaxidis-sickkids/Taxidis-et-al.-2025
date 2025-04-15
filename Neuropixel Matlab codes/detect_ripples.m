function [ripples,rchannel,LFPbp] = detect_ripples(LFP,t)

%--------------------------------------------------------------------------
% Detects ripples based on the RMS of the "Fripple" bandpassed LFP from a
% reference shank/channel, when it goes above a threshold.
% Ripple edges are detected based on a second threshold. Only segments
% corresponding to rat immobility are considered.
% The ripple beginning/ending/peak times are stored as columns in a matrix.
% The ripples are aligned over their peak times (for plotting only).
%
% Created by Jiannis Taxidis, Caltech, USA, February 2013
%--------------------------------------------------------------------------

Fs = 2500;
Fripple = [120,200];                                                        % Ripple frequency range to bandpass LFP

bindur = 26;                                                                 % Bin duration for RMS (10 msec)
[Nch,lt] = size(LFP);
bins = 1 : bindur : (lt-bindur);                                           % Array of locations of bin starting times
time_b = t(bins + bindur/2);
lb = length(bins);

crit1 = 10;
crit2 = 2;

periripdur = 0.1;                                                           % Duration around ripple (300 msec total window length (like Kamran's));
mindur = 0.04;                                                              % Minimum accepted ripple duration (20 msec)
mindist = 0.05;                                                             % Minimum accepted distance between ripples (50 msec)

peririp = ceil(periripdur*Fs);


%% FILTER AND CALCULATE ROOT-MEAN-SQUARE
LFPbp = bandpass(LFP',Fripple,Fs);                                          % Filter the LFP
LFPbp = LFPbp';

RMS = nan(Nch,lb);
for i = 1:lb
    RMS(:,i) = sum(LFPbp(:,bins(i) + [0:bindur]).^2,2);                     % Calculate rms for each bin (IT IS THE SUM OF SQUARES ONLY!!! units = mV^2)
end

%% FIND CHANNEL WITH HIGHEST RMS
[~,rchannel] = max(sum(RMS,2));

RMS = RMS(rchannel,:);
LFPbp = LFPbp(rchannel,:);

%% DETECT ALL SEGMENTS ABOVE THE LFP THRESHOLD AND THEIR LIMITS
stdr = std(RMS);                                                            % Standard deviation of rms
CRIT1 = crit1*stdr;                                                         % Set the detection criterion
CRIT2 = crit2*stdr;                                                         % Set the limits detection criterion
RMS(1) = 0;                                                                 % Set the first and last points to zero
RMS(end) = 0;                                                               % to avoid detecting on the edges.

cross = (RMS' >= CRIT1);                                                     % Find the array locations over the threshold
crossbd = find((cross - circshift(cross,1)) == 1);                          % Find locations whose preceding location is below the threshold (ripple beginning)
crossbd = [crossbd find(cross - circshift(cross,1) == -1)];                 % Add locations whose following location is below the threshold (ripple end)

ripples = zeros(length(crossbd(:,1)),3);
for i = 1:length(crossbd(:,1))                                              % For each ripple beginnning
    rtmp = RMS(crossbd(i,1):crossbd(i,2));                                  % Keep the rms
    rippeak = find(rtmp == max(rtmp),1,'first') + crossbd(i,1);                       % Find the overall location of the rms peak
    ripples(i,1) = find(RMS(1:rippeak) < CRIT2,1,'last');                   % Keep the last location before ripple beginnning below the limit threshold
    ripples(i,2) = find(RMS(rippeak:end) < CRIT2,1,'first') + rippeak - 1;  % and the first location after ripple end below the limit threshold
    ripples(i,3) = rippeak;                                                 % and the location of the ripple peak
end
ripples = time_b(ripples);

[~,r] = unique(ripples(:,1), 'rows');                                       % Keep unique events by checking repeating ripple beginning times
ripples = ripples(r,:);
[~,r] = unique(ripples(:,2), 'rows');                                       % Keep unique events by checking repeating ripple end times
ripples = ripples(r,:);

disp('Initially :')
[~,rdur] = ripple_stats(ripples);

%% REMOVE SHORT AND OVERLAPPING RIPPLES
rshort = (rdur < mindur);                                                   % If any ripple lasts less than 20 msec discard it
ripples(rshort,:) = [];

k = ripples(2:end,1) - ripples(1:(end-1),2);                                % Find the distances between all ripples
k = (k < mindist);                                                          % Find when the distance is not at least 50 msec (covers overlaps too)
k1 = [k; 0];
k2 = [0; k];                                                                % Add a zero for the first ripple that was not included
ripples(k1==1,2) = ripples(k2==1,2);                                        % Unite them (use the ending of the second ripple as end of the first one...
ripples(k2==1,:) = [];                                                      % ...and delete the times of the second one)

disp('After close-ripples removal: ')
ripple_stats(ripples);

%% REMOVE BOUDNARY RIPPLES
boundaries = ripples(:,3)*Fs - peririp <= 0 | ripples(:,3)*Fs + peririp > lt;                      % If the perievent duration hits the LFP boundaries
ripples(boundaries,:) = [];

disp('After boundary-ripples removal: ')
rN = ripple_stats(ripples);

%% PLOT LFPS AND RIPPLE TIMES
% miL = abs(min(LFP(rchannel,:)));
% maLr = max(LFPbp);
% rms = RMS / max(RMS) * maLr;              % Scale RMS to have the same max as the LFPbp for plotting
% z = miL + maLr + max(LFP(rchannel,:));
% 
% figure;
% subplot(211)
% hold on;
% plot(t, LFPbp);
% plot(time_b, rms,'k','LineWidth',1.5);
% line([t(1),t(end)], [CRIT1,CRIT1] / max(RMS) * maLr,'Color','k')
% line([t(1),t(end)], [CRIT2,CRIT2] / max(RMS) * maLr,'Color','k')
% 
% plot(t, LFP(rchannel,:) + miL + maLr,'b');
% 
% for i = 1:rN
%     line([ripples(i,1),ripples(i,1)],[min(LFPbp), z ],'Color','b');
%     line([ripples(i,2),ripples(i,2)],[min(LFPbp), z ],'Color','r');
% end
% axis tight;
%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------



function [rN,rdur] = ripple_stats(ripples)
rN = size(ripples,1);                                                       % Number of ripples in this session
rdur = ripples(:,2)-ripples(:,1);                                           % and ripple durations
rmax = max(rdur);                                                           % Maximum duration
disp([num2str(rN),' ripples of average duration ', num2str(mean(rdur)*1000),' msec and maximum duration ',num2str(rmax*1000),' msec',' (SD = ',num2str(std(rdur)*1000),')'])
