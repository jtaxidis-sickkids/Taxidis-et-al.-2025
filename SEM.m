function S = SEM(A)

% COMPUTE STANDARD ERROR MEAN ACROSS ROWS OF MATRIX A
lc = size(A,2);
S = nan(1,lc);

for c = 1:lc                            % For each column
    N = ~isnan(A(:,c));                 % Find rows without nan entry
    N = sum(N);                         % Count them
    S(c) = nanstd(A(:,c),1)/sqrt(N);    % Compute SEM along column with the appropriate sample size
end

S(isnan(S)| isinf(S)) = 0;