function [asapfile,asapfiles,videonumber] = get_GEVIfile(GEVI,ID,day,session)

drive = pwd;
drive = drive(1);
asapfiles = dir(fullfile([drive,':/',GEVI],ID,day,'*.tif'));                % Find all tif files WITHIN THE SPECIFIC DIRECTORY (not subfolders)

% Remove all motion correction and FOV files
out = [];
for s = 1:length(asapfiles)                                                
    nm = asapfiles(s).name;
    if contains(nm,'_MC') | contains(nm,'FOV') | ...
            contains(nm,'summary_images') | contains(nm,'cell_anatomy') | ...
            contains(nm,'AVG')
       out = [out, s];                                                      
    end
end
asapfiles(out) = [];

asapfiles = struct2cell(asapfiles);
asapfiles = asapfiles(1,:)';

videonumber = zeros(length(asapfiles),1);
for s = 1:length(asapfiles)
	videoname = asapfiles{s};                                               % Get the name of the tif video
    videonumber(s) = str2double(videoname(regexp(videoname,'[0-9]')));      % Find the actual number of the tif video
end

[videonumber,order] = sort(videonumber);                                    % Sort the numbers (avoids e.g. 12 appearing before 2)
asapfiles = asapfiles(order);                                               % Sort videos correctly

for s = 1:length(asapfiles)
    asapfiles{s} =  fullfile([drive,':/',GEVI],ID,day,asapfiles{s});        % Add full directory
end
asapfile = asapfiles{session};
videonumber = videonumber(session);

