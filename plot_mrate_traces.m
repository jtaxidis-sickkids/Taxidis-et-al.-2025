function plot_mrate_traces(M,S,t,tp,cols)

ttypes = size(M,1);
for tt = 1:ttypes
    M(tt,:) = smooth(M(tt,:),5,'moving');                                     % Smooth firing rate
    S(tt,:) = smooth(S(tt,:),5,'moving');                                     % Smooth firing rate
end

onbins = find(t > tp(1) & t < tp(3));

Mon = M(:,onbins) + S(:,onbins);
ymax = max(Mon(:)) + 0.2;

Mon = M(:,onbins) - S(:,onbins);
ymin = min(Mon(:)) - 0.2;

hold on;
shade_timepoints([ymin ymax],tp);                                            % Make shaded rectangles for timepoints
for tt = 1:ttypes       
    fill_plot(t,M(tt,:),S(tt,:),cols(tt,:));                                 % Plot shaded SD rate
end
ylim([ymin ymax]);
xlim([t(1) tp(5)]);
xlabel('time (sec)');
