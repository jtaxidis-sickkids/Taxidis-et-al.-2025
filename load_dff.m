function dff = load_dff(ID,day,session,roi)

asapfile = get_ASAPfile(ID,day,session);
[path,videoname] = fileparts(asapfile);
Datafile = fullfile(path,[videoname,'_data.mat']);
load(Datafile);

dff = V.DFF;
dff(isnan(dff) | isinf(dff)) = 0;
dff(:,1:2,:) = 0;                     % Remove problematic frames
dff = squeeze(dff(roi,:,:));          % Keep only selected ROI (time x trials)
dff = dff';                         % [trials x time]

