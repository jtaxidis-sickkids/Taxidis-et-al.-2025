function plot_FOV(day,lim)

if nargin == 1
    lim = [1,512];
end

[animal,main_folder,day] = get_animal_name(day);
    day = fullfile(main_folder,'Processed_files',day);
load(fullfile(day,'ROIred.mat'));

% C = imfuse(template_g,template_r,'falsecolor','Scaling','independent','ColorChannels',[2 1 0]); % Infuses green template with red template (sets colors as green/red)
C = imfuse(template_g,template_r,'falsecolor','Scaling','independent','ColorChannels','green-magenta'); % Infuses green template with red template (sets colors as green/red)
C = imadjust(C,[0,0.6]);
imshow(C(lim(1):lim(2),lim(1):lim(2),:));                                             % Show only ~350x350µm (360x360 pixels)
hold on;
% line([10, 62],[320 320],'color','w','linewidth',4)                          % Draw line scale of 50µm (52 pixels)