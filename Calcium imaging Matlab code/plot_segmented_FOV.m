function plot_segmented_FOV(day,lim)

if nargin == 1
    lim = [1 512];
end

[~,folder,day] = get_animal_name(day);    
day = fullfile(folder,'Processed_files',day);
load(fullfile(day,'ROIred.mat'));
load(fullfile(day,'ROI.mat'));

F = (template_g(lim(1):lim(2),lim(1):lim(2)));

if max(F(:))/2.2 > min(F(:))    
    imshow(F,[min(F(:)),max(F(:))/2.2]);
else
    imshow(F,[min(F(:)),max(F(:))]);
end
colormap(gray);

CC = CC(keep);                                                              % Keep only the final ROI contours
CC = CC(unmatched_g);
for r = 1:length(CC)
    if all(CC{r} > lim(1)) & all(CC{r} < lim(2))
        plot(CC{r}(1,:)-lim(1),CC{r}(2,:)-lim(1),'Color',[9 117 220]/255, 'linewidth', 1.5); % Plot contours
    end
end

CC_r = [CC_r(matched(:,2)); CC_r(unmatched_r)];
for r = 1:length(CC_r)
    if all(CC_r{r} > lim(1)) & all(CC_r{r} < lim(2))
        plot(CC_r{r}(1,:)-lim(1),CC_r{r}(2,:)-lim(1),'Color','r', 'linewidth', 1.5); % Plot contours
    end
end