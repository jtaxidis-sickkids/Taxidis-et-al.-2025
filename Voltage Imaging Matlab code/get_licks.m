function licks = get_licks(Licks,timepoints)


licktime = (0 : 1e-3 : timepoints(end)-1e-3)';

Ntr = size(Licks,2);
licks = cell(1,Ntr);
Licks = Licks < 1;

for tr = 1:Ntr
    LL = Licks(:,tr) - circshift(Licks(:,tr),1) == 1;
    k = find(LL);
    licks{tr} = licktime(k)';
end