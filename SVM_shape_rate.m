function Rt = SVM_shape_rate(R,cells,bins)

if isrow(bins)
    bins = repmat(bins,length(cells),1);
end

Rt = zeros(length(cells),size(bins,2),size(R,3));
for c = 1:length(cells)
    Rt(c,:,:) = R(cells(c),bins(c,:),:);                                    % Keep only selected cells and bins
end
Rt = squeeze(mean(Rt,2));                                                   % Mean Rate over the bins (cells x trials);
if length(cells) == 1; Rt = Rt';end                                         % If there is one cell, this will fix the squeeze-transposition
Rt = zscore(Rt,[],2);                                                       % Z-score over trials
Rt = Rt';                                                                   % trials x cells