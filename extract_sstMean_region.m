% EXTRACT MEAN SST IN PREDEFINED REGION FOR JOST'S ANALYSIS
clear

% Data files
years = 2006:2015;
basedatapath = '/Volumes/rap/martigalitapias/ERA-Interim/SST_46km/'; % on taku-leifr
grid = dlmread('~/Desktop/Grids_maps/grids/grid46km_45N.txt');
lat = grid(:,2);
lon = grid(:,1);

% Create index list to extract sst in region of interest.
% Looking at Jost's figures, the region expands down to 55N in the Atlantic
% sector, but only down to 65N in the remaining sectors.
indlist = find(lat > 55 & ((lon > -60 & lon < 40) | (lat > 65 & (lon < -60 | lon > 40))));

% Select latitudes within region
latin = lat(indlist);

tic
for yyyy = years
    
    datapath = [basedatapath num2str(yyyy) '/']; % uncomment for taku-leifr
    flist = dir([datapath '/*.txt']);
    
    for mmdd = 1:length(flist)
        
        disp(['Opening ' flist(mmdd).name])
        filepath = ([datapath flist(mmdd).name]);
        datain = dlmread(filepath);
        datain = datain(indlist);
        datain(datain==-999) = nan;
        
        % Write to summary file
        dummy1 = regexp(flist(mmdd).name,'\.txt','split');
        date_char = char(dummy2{1,1});
        
        tmp = [yyyy...
            str2double(date_char(9:10))...
            str2double(date_char(11:12))...
            nanmean(datain)... % mean sea surface temperature in selected region
            nanmin(datain)... % min sea surface temperature in selected region
            nanmax(datain)... % max sea surface temperature in selected region
            sum(~isnan(datain))... % number of 46 km pixels used to calculate the mean regional sst
            mean(latin(~isnan(datain)))... % mean latitude of 46 km pixels used to calculate the mean regional sst
            max(latin(~isnan(datain)))]; % max latitude of 46 km pixels used to calculate the mean regional sst
        
        dlmwrite('./summary_sst_SpitsbergenRegion.txt',tmp,'-append');
        
    end
end

toc
