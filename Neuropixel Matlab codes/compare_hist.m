function compare_hist(X,C1,C2)
X1 = X(C1);
X2 = X(C2);
m1 = min(X)-0.1;
m2 = max(X)+0.1;
binstep = (m2-m1)/25;                                                       % 30 bins
histogram(X1,m1:binstep:m2,'Facecolor','k');   
histogram(X2,m1:binstep:m2,'Facecolor','r');   
axis tight;