function Plot_av_spectrogram(INclass)

%% LOAD ALL SPECTROGRAMS
if strcmp(INclass,'PV')
    Days_PV;
    maxd = 7;
elseif strcmp(INclass,'SOM')
    Days_SOM;
    maxd = 8;
end

Na = length(asap);
PSD = cell(Na,maxd);

for a = 1:Na
    ID = asap(a).name;
    days = asap(a).sessions(:,1);
    ses = asap(a).sessions(:,3);
    for d = 1:length(days)                                               % For each selected day      
        PSD{a,d} = [];
        
        for s = 1:length(ses{d})                                    % for each selected session that day
            [asapfile,~] = get_ASAPfile(ID,days{d},ses{d}(s));
            [path,videoname] = fileparts(asapfile);
            disp(asapfile);
 
            datafile = fullfile(path,[videoname,'_UnitAnalysis.mat']);
            load(datafile,'F');           

            PSD{a,d} = cat(1,PSD{a,d},{F.PSD});                       
        end
    end
end

fq = F.fq;
tq = F.tq;

%% STACK IN SINGLE COLUMNS & BREAK UP THE DATA CONTAINING MULTIPLE CELLS
PSD = vertcat(PSD{:});

count = 1;
PSDall = [];
for i = 1:length(PSD)
    for r = 1:size(PSD{i},1)
        PSDall{count,1} = squeeze(PSD{i}(r,:,:,:));                          % Keep only that cell's rates [time x freq x trials]
        count = count + 1;
    end
end
PSD = PSDall;

clear PSDall 

%% GET MEAN PSD ACROSS TRIALS PER CELL
mPSD = cellfun(@(x) nanmean(x,3), PSD, 'UniformOutput',false);
mPSD = reshape(mPSD,[1,1,length(mPSD)]);
mPSD = cell2mat(mPSD);
mPSD = 10*log10(mPSD); % Turn to decibels)M
   
% REMOVE 1/f COMPONENT TO FLATTEN IT
mmPSD = nanmean(mPSD,3);

k = mean(mmPSD,2);
fout = fit(fq'+0.001,k,'power2'); 
y = fout.a*((fq+0.001).^fout.b)+fout.c;
y = repmat(y',[1,length(tq)]);
mmPSD = mmPSD - y;

%% PLOT
figure;

% KEEP UP TO 18Hz  
fi = find(fq <= 18,1,'last');
fq(fi+1:end) = [];

k = mmPSD(1:fi,:);
% ----

imagesc(tq,fq(1:fi),k,[-4.5 4.5]);
colorbar('off');
colormap parula
set(gca,'YDir','normal');
line([1 1],[0 fq(fi)],'Color','w');
line([2 2],[0 fq(fi)],'Color','w');
line([7 7],[0 fq(fi)],'Color','w');
line([8 8],[0 fq(fi)],'Color','w');
xlim([0.7 10.3]);



