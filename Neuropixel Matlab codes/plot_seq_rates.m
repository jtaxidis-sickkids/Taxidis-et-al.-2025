function plot_seq_rates(M,cells,btime,tp,maxcol)

hold on;

%% SMOOTH
[Ns,Nt] = size(M);                                                          % Cells x bins

% M = smoothdata(M,2,'gaussian',12);                                    % Smooth its firing rate

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

imagesc(btime , 1:Ns, M, [0 maxcol]);                                        % Plot bins x cells x rate
set(gca,'YDir','reverse');                                              % START PLOTTING FROM THE TOP
for i = 1:length(tp)-1
    line(tp(i)*[1 1],[0.5 Ns+0.5], 'Color','w','linewidth',1.5);        % Make lines for timepoints
end
set(gca,'Ytick',ystep:ystep:Ns,'YTickLabel',{ystep:ystep:Ns});
ylim([0.5 Ns+0.6]);


% if length(tp) >= 7
%     xlim([btime(1) tp(7)]);
% else
%     xlim([btime(1) tp(3)]);
% end
xlabel('time (sec)');
ylabel('Mean normalized rates');
axis tight
