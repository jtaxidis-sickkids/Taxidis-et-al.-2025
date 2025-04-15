function p = powerlaw_goodness_of_fit(x,y,f)

reps = 500;

lx = length(x);

rng(0);
rn = rand(reps,lx);

%% GET DISTRIBUTION OF FITTED DATA AND KOLMOGOROV-SMIRNOV STATISTIC
yfit = f.a * x.^f.b + f.c;
KS = max(abs(cumsum(y) - cumsum(yfit)));

%% COMPUTE CDF OF NORMALIZED PROBABILITY DISTRIBUTION
pmf = yfit / sum(yfit);
cdf = cumsum(pmf);

%% DRAW RANDOM POINTS FROM THE SAME PROBABILITY DISTRIBUTION, FIT THEM, AND GET THEIR KS STATISTIC
KSr = zeros(1,reps);
xr = nan(1,lx);
for r = 1:reps 
    for i = 1:lx                                                            % For each of the original points
        if rn(r,i) == 1                             
            xr(i) = max(x);
        else
            idx = sum(cdf < rn(r,i)) + 1;                                   % Find the index of the random number in the distribution
            xr(i) = x(idx);                                                 % Keep the corresponding point of the distribution
        end
    end 
    yr = hist(xr,x);                                                        % Get their distribution
    yr = 100 * yr/lx;                                                       % Turn to % ratio

    fr = fit(x',yr','power2');                                              % Fit with power law
    yrfit = fr.a * x.^fr.b + fr.c;                                          % Get data of fit
    KSr(r) = max(abs(cumsum(yr) - cumsum(yrfit)));                          % and compute KS statistic
end
hold off;

p = sum(KSr < KS) / reps;                                                   % Estimate p-value
