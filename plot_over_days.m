function plot_over_days(A)

hold on;

G = lines(1);
lineC = [80 80 80]/255;

[Na,ld] = size(A);                                                             % Number of animals and days
dayaxis = 1:ld;

markers = ['o';'s';'d';'^';'v';'*';'+';'x';'>';'<';'.';'o';'o';'o';'o']; 
if Na > length(markers)
    markers = repmat('.',[Na 1]);
end

mA = nanmean(A,1);                                                          % Get mean performance and SE
sA = SEM(A);

for a = 1:Na                                                                % For each animal
    plot(dayaxis,A(a,:),'-','Marker',markers(a),'Color',lineC,'MarkerSize',3);                                            % Plot Performance per day for Step2
end
fill_plot(dayaxis,mA,sA,G);
axis tight;