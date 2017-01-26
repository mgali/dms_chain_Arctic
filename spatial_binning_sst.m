% SPATIAL BINNING OF ERA-Interim SST DATA
clear

dateoutformat = 'YYYYDDD'; % YYYYDDD or YYYYMMDD
outformat = 'txt'; % 'nc' or 'txt'

% Grids
kmgrid2 = '28'; % 28, 37 or 46
gridpath = '~/Desktop/Grids_maps/grids/'; % TESTS and taku-leifr
grid1name = 'ERAi_45N.mat';
load([gridpath grid1name]); grid1 = [lon lat]; clear lon; clear lat
grid2name = ['grid' kmgrid2 'km_45N.txt'];
grid2 = dlmread([gridpath grid2name],' ');
npixels2 = size(grid2,1);

% Data files
years = 1998:2015; % 1998:2015;
varnames = {'sst'};
% basedatapath = '~/Desktop/Geophysical_data/ERA-Interim/SST/'; baseoutpath = basedatapath; % FOR TESTS
basedatapath = '/Volumes/rap/martigalitapias/ERA-Interim/SST/'; % on taku-leifr
baseoutpath = sprintf('/Volumes/rap/martigalitapias/ERA-Interim/SST_%skm/',kmgrid2); % on taku-leifr

%% Use index list to bin dmspt from grid1 to grid2

load(['indlist_ERAi_45N_to_' kmgrid2 'km.mat'])
tic
nrows = size(indlist,1);
ncols = length(varnames); % usual definition
% ncols = length(days); dataout = nan(nrows,ncols); % TESTS, alternative to compare consecutive images

toc
for vv = 1:length(varnames)
    for yyyy = years
        
        datapath = [basedatapath num2str(yyyy) '/'];
        flist = dir([datapath '/*.nc']);
        
        % Select days (eg first of 8D period)
        % flist(mod(1:size(flist,1),8)~=1) = [];
        
        for mmdd = 1:length(flist)
            
            disp(['Opening ' flist(mmdd).name])
            filepath = ([datapath flist(mmdd).name]);
            datain = ncread(filepath,varnames{vv});
            datain = datain(:) - 273.15;
            
            % Binning/reprojection
            dataout = repbin_grid1_grid2(datain,indlist);
            
            % Write file
            dataout(isnan(dataout)) = -999;
            dummy1 = regexp(flist(mmdd).name,'\.nc','split');
            dummy2 = regexp(dummy1{1,1},'\ERAi_','split');
            datename = dummy2{1,2};
            
            if strcmp(dateoutformat,'YYYYMMDD')
                outpath = sprintf('%s%i/%s_%s.%s',baseoutpath,yyyy,varnames{vv},datename,outformat); % taku-leifr
            elseif strcmp(dateoutformat,'YYYYDDD')
                Y = str2double(datename(1:4));
                M = str2double(datename(5:6));
                D = str2double(datename(7:8));
                DOY = yearday(D,M,Y,0,0,0);
                outpath = sprintf('%s%0.0f/%s_%0.0f%03.0f.%s',baseoutpath,Y,varnames{vv},Y,DOY,outformat); % taku-leifr
            end
            
            if strcmp(outformat,'txt')
                dlmwrite(outpath,dataout,'\t');
            elseif strcmp(outformat,'nc')
                nccreate(outpath,varnames{vv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                ncwrite(outpath,varnames{vv},dataout);
            end
            
            % Store mean temperature in separate file
            % tmp = [yyyy mmdd nanmean(dataout)];
            % dlmwrite('./summary_MeanSstDay.txt',tmp,'-append');
            
            % % Plot side by side datain and dataout, save
            % h = figure(); clf, set(gcf,'units','centimeters','position',[5 20 40 20])
            % subplot(1,2,1)
            % plot_sst_stereo(datain,8,grid1)
            % subplot(1,2,2)
            % % dataout(isnan(dataout)) = 10;
            % plot_sst_stereo(dataout,8,grid2)
            % % print(h,sprintf('~/Desktop/Heintzenberg_Leck_collaboration/SST_maps/SST_%s.png',datename),'-dpng');
            % print(h,sprintf('%sSST_%s_ERAvs%skm.png',datapath,datename,kmgrid2),'-dpng');
            % close(h)
            
        end
    end
end
toc
