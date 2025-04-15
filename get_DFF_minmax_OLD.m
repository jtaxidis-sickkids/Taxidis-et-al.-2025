function [MinD,MaxD,oddff] = get_DFF_minmax(dff,time,timepoints)
    
basedff = dff(:,time >= 0.5 & time <= timepoints(1));                   % BASELINE
mu = std(basedff,[],2);                                                 % STD ACROSS TIME FOR EACH TRIAL
sigma = std(basedff,[],2);                                               % STD ACROSS TIME FOR EACH TRIAL

oddff = dff(:,time >= timepoints(1) & time <= timepoints(2));           % KEEP ONLY ODOR BINS
oddff = (oddff - mu) ./ sigma;                                          % Z-SCORE

MinD = min(oddff,[],2);                                          % Find minimum over odors for each trial
MaxD = max(oddff,[],2);                                          % Find minimum over odors for each trial
