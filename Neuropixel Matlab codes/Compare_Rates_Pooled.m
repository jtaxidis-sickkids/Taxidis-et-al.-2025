function Compare_Rates_Pooled(MR,bins)

titls = {'PY';'IN'};
cols = ['k';'b';'r'];
       
od1bins = find(bins > 1 & bins < 2);
lbod = length(od1bins);

MRod = cell(size(MR));
for i = 1:4
    MRod{i} = nanmean(MR{i}(:,od1bins),2);                                  % Mean zscored signal over the odor
end
Nr = cellfun(@length,MRod(1,:));

%% COMPARE POOLED RATES WITHvsWITHOUT OPTO
figure;
for ct = 1:2
    subplot(2,5,(ct-1)*5+1); hold on;
    
    k = MRod{1,ct} >= 1;                             % Index cells with >0 odor response

    % PLOT MEAN ZSCORED SIGNALS ACROSS ALL CELLS OR WITH POSITIVE/NEGATIVE ODOR RESPONSE AT LIGHT OFF
    fill_plot(bins,nanmean(MR{1,ct},1),SEM(MR{1,ct}),'k');
    fill_plot(bins,nanmean(MR{1,ct}(k,:),1),SEM(MR{1,ct}(k,:)),'b');
    fill_plot(bins,nanmean(MR{1,ct}(~k,:),1),SEM(MR{1,ct}(~k,:)),'r');
    
    % REPEAT FOR LIGHT ON TRIALS WITH DASHED LINES
    fill_plot(bins,nanmean(MR{2,ct},1),SEM(MR{2,ct}),'k'); 
    fill_plot(bins,nanmean(MR{2,ct}(k,:),1),SEM(MR{2,ct}(k,:)),'b');
    fill_plot(bins,nanmean(MR{2,ct}(~k,:),1),SEM(MR{2,ct}(~k,:)),'r');
    hline = findobj(gcf, 'type', 'line');set(hline(1:3),'LineStyle','--');  % Set lines as dashed

    axis tight;
    xlim([0.5,2.5]);            % PLot only around odor
    title(titls{ct});

    for j = 1:3
        ktemp = ones(Nr(ct),1) * (j == 1) + ...
            k * (j == 2) + ...
            ~k * (j == 3);
        ktemp = logical(ktemp);

        M1 = MR{1,ct}(ktemp,od1bins);
        M2 = MR{2,ct}(ktemp,od1bins);
     
        if sum(ktemp)>0
            pv = nan(1,lbod);
            for b = 1:lbod
               pv(b) = ranksum(M1(:,b), M2(:,b));
            end
            [~,pv] = fdr(pv);
            
            plot(bins(od1bins(pv < 0.05)), (3+0.2*j)*ones(1,sum(pv < 0.05)),'o','Color',cols(j));
        end
    end
end

%% COMPARE MEAN ODOR RESPONSES
for ct = 1:2
    k = MRod{1,ct} >= 1;                             % Index cells with >0 odor response
    for j = 1:3
        subplot(2,5,(ct-1)*5 + 1+j); hold on;

        ktemp = ones(Nr(ct),1) * (j == 1) + ...
            k * (j == 2) + ...
            ~k * (j == 3);
        ktemp = logical(ktemp);
        
        M1 = MRod{1,ct}(ktemp);
        M2 = MRod{2,ct}(ktemp);
     
        if ~isempty(M1)
        plot([1 2],[M1,M2],'o-','Color',cols(j));
        errorbar([0.9 2.1],[nanmean(M1) nanmean(M2)],[SEM(M1) SEM(M2)],'o-k')
        set(gca,'Xtick',[1,2],'XTickLabel',{'no opto';'opto'});
        ylabel('Mean odor responses');
        
        [~,pM] = ttest(M1,M2)
        plot_significance(pM,0.9,2.1,nanmean(M1), nanmean(M2))
        end
    end
end

%% COMPARE DISTRIBUTIONS OF ODOR RESPONSES WITHvsWITHOUT OPTO
for ct = 1:2
    subplot(2,5,ct*5); hold on;
       
    binstep = 4*(ct==1) + 2*(ct==2);
    edges = 1 : binstep : Nr(ct);                                           % Make bins
    de = edges(2)-edges(1);                                             % bin size
    
    for i = 1:2
        plot(1:Nr(ct),MRod{i,ct},'.','Markersize',4);
        line([1 Nr(ct)],[1 1],'LineStyle','--','Color','k');
        
        [M,S] = bin_x_axis(1:Nr(ct),MRod{i,ct},edges);                          % Bin the MEAN ZSCORED SIGNAL
        fill_plot(edges(1:end-1)+de/2, M, S, 'b');                          % Plot
    end
    hline = findobj(gca,'type','line');
    set(hline(1),'LineStyle','--');
    axis tight;
    
    [~,pks] = kstest2(MRod{1,ct},MRod{2,ct});
%     pks = signrank(MRod{1,ct},MRod{2,ct});
    title(num2str(pks));
    %
    %     for i = 1:2
    %         mdl = fitlm(1:Nr,MRod{i,ct})
    %         tbl = anova(mdl)
    %     end
    %
    for i = 1:2
        [ffit,~,out] = fit((1:Nr(ct))',MRod{i,ct},'power2');
        ffit
        plot(ffit,(1:Nr(ct))',MRod{i,ct});
        
        %         p = powerlaw_goodness_of_fit((1:Nr),MRod{i,ct}',ffit)
        
    end
end


