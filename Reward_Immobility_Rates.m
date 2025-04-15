function Ri_rew = Reward_Immobility_Rates(INclass)

%% INITIALIZE
if strcmp(INclass,'PV')
    Days_PV;
elseif strcmp(INclass,'SOM')
    Days_SOM;
end

%% INITIALIZE
Na = length(asap);

Ri_rew = [];

%% POOL FILES
for a = 1:Na                                                                % For each selected animal
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1:length(days)                                               % For each selected day
        for s = 1:length(ses{d})                                    % for each selected session that day
            asapfile = get_ASAPfile(ID,days{d},ses{d}(s));
            [path,videoname] = fileparts(asapfile);
            datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
            load(datafile,'M','V');                                 % Load the Theta-Motion analysis file
            disp(asapfile);
            
            bins = V.bins;
            R = V.R;
            Tm = M.Move_time;
            clear M V
            
            lb = length(bins);
            
            bin_rew = find(bins < 9 , 1, 'last');
            
            Riw = [];
            for tr = 1:size(R,3)
                bi = true([1,lb]);
                bi(1:bin_rew) = 0;
                
                if ~isempty(Tm{tr})
                    mot = [];
                    for m = 1:size(Tm{tr},1)                                % For each motion segment
                        [~,k1] = min(abs(bins - Tm{tr}(m,1)));              % Find the bin closest to motion initiation
                        [~,k2] = min(abs(bins - Tm{tr}(m,2)));              % and termination
                        k2 = min(k2,lb);                                    % In case it goes just above lb
                        mot = [mot, k1:k2];                                 % Keep the motion indexes
                    end
                    bi(mot) = 0;
                end
                
                Riw = cat(2,Riw,R(:,bi,tr));
            end
            Ri_rew = [Ri_rew; nanmean(Riw,2)];
        end
    end
end
