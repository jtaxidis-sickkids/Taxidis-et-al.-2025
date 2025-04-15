function plot_seq_rates(M,cells,btime,tp,Cmax)

hold on;

%% SMOOTH 
if ~isempty(cells)
    M = M(cells(:,1),:);                                                        % Keep the average rate of the M-cells
end
[Ns,Nt] = size(M);                                                          % Cells x bins

for c = 1:Ns                                                                % For each M-cell
    M(c,:) = smooth(M(c,:),5,'moving');                                     % Smooth its firing rate
end

lb = length(btime);
if Nt < lb
    M = [M, zeros(Ns,lb-Nt)];                                               % If missing time points, add zeros
else
    M = M(:,1:lb);                                                          % Keep only the btime entroes (in some cases M has more/less entries)
end

%% PLOT
if Ns < 20
    ystep = 1;
elseif Ns >=20 & Ns < 50
    ystep = 5;
elseif Ns >= 50 & Ns < 300
    ystep = 20;
elseif Ns >= 300 & Ns < 500
    ystep = 50;
elseif Ns >= 500
    ystep = 100;
end

% color-plot
imagesc(btime , 1:Ns, M, [-Cmax/2 Cmax]);                                        % Plot bins x cells x rate
set(gca,'YDir','reverse');                                              % START PLOTTING FROM THE TOP
for i = 1:length(tp)-1
    line(tp(i)*[1 1],[0.5 Ns+0.5], 'Color','w','linewidth',1.5);        % Make lines for timepoints
end
set(gca,'Ytick',ystep:ystep:Ns,'YTickLabel',{ystep:ystep:Ns});
ylim([0.5 Ns+0.6]);
if length(tp) >= 7
    xlim([btime(1) tp(6)]);
else 
    xlim([btime(1) tp(3)]);
end
% xlabel('time (sec)');
% ylabel('Mean normalized rates');
