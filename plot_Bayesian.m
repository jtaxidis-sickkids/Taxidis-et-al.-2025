function [T,TT,B] = plot_Bayesian(T,TT,B,Tsh,TTsh,Bsh,bins,sub1,sub2)
    
%% SHAPE DATA
Bt = B;
for i = 1:length(T)                                     % For each decoded cell                          
    lpred = size(T{i},1);                               % Count the number of decoded trials
    Bt{i} = repmat(Bt{i},lpred,1);                      % Repeat the decoded-bins vector that many times (to have same vector sizes)
    Bt{i} = Bt{i}(:)';                                   % Stack them all in one row
    T{i} = T{i}(:)';
%     Bt{i} = Bt{i}';
%     T{i} = T{i}';
end

Btsh = Bsh;
for i = 1:length(Tsh)                                   % Repeat for shuffles (same process x reps)
    lpred = size(Tsh{i},1);
    Btsh{i} = repmat(Btsh{i},lpred,1);
    Btsh{i} = Btsh{i}(:)';
    Tsh{i} = Tsh{i}(:)';
%     Btsh{i} = Btsh{i}';
%     Tsh{i} = Tsh{i}';
end


B = cell2mat(B);                                                         % Keep decoded bins of all days
Bt = cell2mat(Bt);                                                         % Keep decoded bins of all days
T = cell2mat(T);                                                            % and time prediction errors
TT = cell2mat(TT);                                                          % and trial prediction accuracy (can contain Nans for non-decoded bins in trials)
Bsh = cell2mat(Bsh);                                                     %
Btsh = cell2mat(Btsh);                                                     %
Tsh = cell2mat(Tsh);                                                        % Same for shuffled-cells trials
TTsh = cell2mat(TTsh);                                                      %


T = abs(T);                                                                 % Keep absolute values of time-prediction error
Tsh = abs(Tsh);

lb = length(bins);

%% PLOT DECODING TIME ERRORS ----------------------
subplot(sub1); hold on;
plot(Bt+rand(size(Bt))*0.1 , T+rand(size(T))*0.1 ,'.','Color',[.8 .8 .8]);    % Plot ALL trial prediction accuracies (plus noise)
[M,S,D] = bin_x_axis(Bt,T,bins);
fill_plot(bins(1:end-1),M,S,'k');

[Mb,~,Db] = bin_x_axis(Btsh,Tsh,bins);
plot(bins(1:end-1),Mb,'--r');                                                        % Plot shuffled cells baseline

xlabel('time (sec)');ylabel('Time decoding error (sec)');
xlim([1 7]); ylim([0 max(T(:))]);

%% PLOT SIGNIFICANCE
p = nan(1,lb);
for b = 1:lb-1                                                                % For each histogram bin
    p(b) = significance(D{b},Db{b},'smaller');                                % Test significance of trial prediction over shuffled prediction
end
[~,p] = fdr(p);                                                             % FDR Correction
plot(bins(p < 0.05),ones(1,sum(p < 0.05))*max([M,Mb])*1.1,'*k');

%% PLOT DECODING TRIAL TYPE ACCURACIES ----------------------
subplot(sub2); hold on;
% plot_Bayesian(Btt,TT,Bttsh,TTsh,bins,'larger');                              % Plot decoded trial type accuracy and significance
plot(B+rand(size(B))*0.1 , TT+rand(size(TT))*0.1 ,'.','Color',[.8 .8 .8]);    % Plot ALL trial prediction accuracies (plus noise)
[M,S,D] = bin_x_axis(B,TT,bins);
fill_plot(bins(1:end-1),M,S,'k');

[Mb,~,Db] = bin_x_axis(Bsh,TTsh,bins);
plot(bins(1:end-1),Mb,'--r');                                                        % Plot shuffled cells baseline

xlabel('time (sec)');ylabel('Accuracy of trial-type decoding per session (%)');
xlim([1 7]); ylim([40 110]);

%% PLOT SIGNIFICANCE
p = nan(1,lb);
for b = 1:lb-1                                                                % For each histogram bin
    p(b) = significance(D{b},Db{b},'larger');                                % Test significance of trial prediction over shuffled prediction
end
[~,p] = fdr(p);                                                             % FDR Correction
plot(bins(p < 0.05),ones(1,sum(p < 0.05))*max([M,Mb])*1.1,'*k');


