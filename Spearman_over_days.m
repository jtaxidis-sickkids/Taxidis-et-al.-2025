function [R,P] = Spearman_over_days(A)


[Na,days] = size(A);

D = repmat(1:days,Na,1);

D = D(:);
A = A(:);

D(isnan(A)) = [];
A(isnan(A)) = [];

[R,P] = corr(D,A,'type','Spearman');

