function mark_field(pr_od,mbin,trials)    

Ntr = length(trials);

if pr_od == 1
    line([1 1]*mbin,[1 sum(trials == 1)],'color','w','linestyle','--','linewidth',2); % Add dashed line on its field bin
elseif pr_od == 2
    line([1 1]*mbin,[sum(trials == 1)+1 Ntr],'color','w','linestyle','--','linewidth',2); % Add dashed line on its field bin
elseif pr_od == 0
    line([1 1]*mbin,[1 Ntr],'color','w','linestyle','--','linewidth',2); % Add dashed line on its field bin
end