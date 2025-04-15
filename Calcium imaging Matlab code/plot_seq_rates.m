function plot_seq_rates(M,S,cells,onbins,btime,tp,flag)

hold on;

%% SMOOTH 
if ~isempty(cells)
    M = M(cells(:,1),:);                                                        % Keep the average rate of the M-cells
end
[Ns,Nt] = size(M);                                                          % Cells x bins

for c = 1:Ns                                                                % For each M-cell
    M(c,:) = smooth(M(c,:),5,'moving');                                     % Smooth its firing rate
end

if flag == 1                                                                % If plotting traces
    if ~isempty(cells)
        S = S(cells(:,1),:);                                                    % Keep the SD rate
    end
    for c = 1:Ns                                                                % For each M-cell
        S(c,:) = smooth(S(c,:),5,'moving');                                     % Smooth its firing rate
    end
    
    if Ns < 5                                                               % Also set traces colors
        C = cell_colors;
    else
        C = parula(Ns);
    end
end

lb = length(btime);
if Nt < lb
    M = [M, zeros(Ns,lb-Nt)];                                               % If missing time points, add zeros
    S = [S, zeros(Ns,lb-Nt)];
else
    M = M(:,1:lb);                                                          % Keep only the btime entroes (in some cases M has more/less entries)
    if flag == 1
        S = S(:,1:lb);
    end
end

%% PLOT
if Ns < 20
    ystep = 1;
elseif Ns >=20 & Ns < 50
    ystep = 5;
elseif Ns >= 50 & Ns < 300
    ystep = 20;
elseif Ns >= 300 & Ns < 1000
    ystep = 100;
elseif Ns >= 1000 & Ns < 10000
    ystep = 200;
elseif Ns >= 10000
    ystep = 1000;
end

if flag == 1                                                                % If traces-plot
    shade_timepoints([0 Ns+2],tp);                                          % Make shaded rectangles for timepoints
    for c = 1:Ns                                                            % For each cell
        fill_plot(btime,M(c,:) + (Ns-c+1),S(c,:),'k');                    % Plot shaded SD rate
        plot(btime,M(c,:) + (Ns-c+1),'Linewidth',1.5,'Color',C(c,:));       % And mean Rate (STARTING FROM THE TOP OF FIGURE)
    end
    set(gca,'Ytick',ystep:ystep:Ns,'YTickLabel',{fliplr(ystep:ystep:Ns)});
    ylim([0.5 Ns+2]);
    
elseif flag == 2                                                            % If color-plot
    imagesc(btime , 1:Ns, M, [0 1]);                                        % Plot bins x cells x rate
    set(gca,'YDir','reverse');                                              % START PLOTTING FROM THE TOP
    for i = 1:length(tp)-1
        line(tp(i)*[1 1],[0.5 Ns+0.5], 'Color','w','linewidth',1);        % Make lines for timepoints
    end
    set(gca,'Ytick',ystep:ystep:Ns,'YTickLabel',{ystep:ystep:Ns});
    ylim([0.5 Ns+0.6]);
end

if length(tp) >= 7
    xlim([btime(1) tp(7)]);
else 
    xlim([btime(1) tp(3)]);
end
% xlabel('time (sec)');
% ylabel('Mean normalized rates');
