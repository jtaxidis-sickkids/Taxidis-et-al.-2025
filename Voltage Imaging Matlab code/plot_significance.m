function plot_significance(p,x1,x2,A1,A2)

A = [A1(:);A2(:)];
A = max(A);
incr = A*0.2;

ml = mean([x1,x2]);

if p < 0.05
    
    line([x1,x1] , A+[0, incr],'Color','k');
    line([x2,x2] , A+[0, incr],'Color','k');
    line([x1,x2] , (A+incr)*[1 1],'Color','k');
    
    if p > 0.01
        plot(ml , A+incr*2,'*k','MarkerSize',5);
    elseif p <= 0.01 && p > 0.001
        plot(ml + [-ml/50,ml/50] , (A+incr*2)*[1 1] , '*k','MarkerSize',5);
    elseif p <= 0.001
        plot(ml + [-ml/25,0,ml/25] , (A+incr*2)*[1 1 1] , '*k','MarkerSize',5);
    end
    
end