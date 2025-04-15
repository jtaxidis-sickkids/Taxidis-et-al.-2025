function [earlylicks,licks] = get_early_lick_trials(day,timepoints)

[~,folder,dd] = get_animal_name(day);
load(fullfile(folder,'Processed_files',dd,'Behavior.mat'));

licktime = [0 : 1e-3 : timepoints(end)-1e-3]';
earlylicks = [];

nomat = cellfun(@isempty, Licks);                                           % FOR ONE CASE OF A MISSING MATLAB FILE (JT77 3/6/17)
Licks(nomat) = [];

licks = cell(length(Licks),20);

for s = 1:length(Licks)
    L = Licks{s};
    L = L < 1;
    
    for tr = 1:size(L,2)
        LL = L(:,tr) - circshift(L(:,tr),1) == 1;
        k = find(LL);
        licks{s,tr} = licktime(k)';
        early = sum(licks{s,tr} >= timepoints(3) & licks{s,tr} < timepoints(4));
        if early > 0
            earlylicks = [earlylicks (s-1)*20 + tr];
        end
    end
end

licks = licks';
licks = licks(:);                       % Stack licktimes by trial
