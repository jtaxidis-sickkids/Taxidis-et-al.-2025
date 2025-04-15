function [M,S,D] = bin_x_axis(X,Y,edges)

lb = length(edges)-1;
M = zeros(1,lb);
S = zeros(1,lb);
D = cell(1,lb);

[Hcd,~,B] = histcounts(X,edges);                        % Count X points within each bin and get the bin where each datapoint belongs

for b = 1:lb                                            % For each bin
    M(b) = nanmean(Y(B == b));                          % Get the mean Y within that spatial bin
    S(b) = nanstd(Y(B == b)) / sqrt(Hcd(b));            % and its SE
    D{b} = Y(B == b);                                   % and all the Y values in that bin
end