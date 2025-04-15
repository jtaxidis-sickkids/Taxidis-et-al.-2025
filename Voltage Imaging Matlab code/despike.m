function [DespDFF,params] = despike(DFF,spiketimes,Fs)

Nt = length(DFF);

%Find a good value for g according to the method shown in the appendix,
%fitting to a function with a modest number of free parameters
g = fitLFPpowerSpectrum(DFF,0.5,250,Fs);

% Assume that a spike lasts from -5 ms to +19 samples (25 ms total) compared to its peak
pre = 5;
post = 19;
Bs = eye(pre+post+1);

% Make a spike entries vector
spiketimes = spiketimes*Fs + 1;                                             % Turn to datapoints
spiketimes = round(spiketimes);
spiketimes = spiketimes + 2;                                                % Add 2 because the actual peak is 2 point after
spiketimes(spiketimes-pre < 1) = [];

S = zeros(Nt,1);
S(spiketimes - pre) = 1;

% Finally despike the LFP
opts.displaylevel = 0;
results = despikeLFP(DFF,S,Bs,g,opts);

% Arrange output
DespDFF = results.z;
params = rmfield(results,'z');