function [edrfile,matfile] = get_Behaviorfiles(ID,day,session,asapfile,asapfiles)

drive = pwd;
drive = drive(1);

edrfiles = subdir(fullfile([drive,':/Jiannis/ASAP/'],ID,day,'*.EDR'));                          % Find all raw files (full pathname)
edrfile = edrfiles(session).name;

matfiles = subdir(fullfile([drive,':/Jiannis/ASAP'],ID,day,'*.mat'));                          % Find all raw files (full pathname)
out = [];
for s = 1:length(matfiles)
    nm = matfiles(s).name;
%     if any([contains(nm,'_SH'), contains(nm,'_ROI'), contains(nm,'_DFF'),...
%             contains(nm,'_Spikes'), contains(nm,'_Behavior'), contains(nm,'_data'),contains(nm,'_Mcells')]) 
    if any([contains(nm,'_data'), contains(nm,'_UnitAnalysis'),contains(nm,'_Mcells')]) 
        out = [out, s];                                                      % Store its number
    end
end
matfiles(out) = [];

% DOESNT WORK BECAUSE FILES DO NOT CARRY CORRECT CREATION TIMESTAMP
% for s = 1:length(matfiles)
%     matfiletimes(s) = datenum(matfiles(s).date);
% end
% --------------------------
matfiletime = nan(1,length(matfiles));
for s = 1:length(matfiles)
    matfilename = matfiles(s).name;
    k = strfind(matfilename,'_');
    k = k(end-1);
    mattime = matfilename(k+1 : end-4);
    
    k = strfind(mattime,'_');
    hours = str2num(mattime(1:k-1));
    minutes = str2num(mattime(k+1:end));
    matfiletime(s) = hours*60+minutes;
end
[~,order] = sort(matfiletime);
matfiles = matfiles(order);
matfile = matfiles(session).name;

[~,nm] = fileparts(asapfile);
disp(['ASAP file: ',nm])
[~,nm] = fileparts(edrfile);
disp(['EDR file: ',nm])
[~,nm] = fileparts(matfile);
disp(['MAT file: ',nm])
if length(asapfiles) ~= length(matfiles) || length(edrfiles) ~= length(matfiles)                                                % Number of recording sessions
    warning('ASAP, EDR or Mat files missing!');
end

