% SPATIAL BINNING TESTS WITH MATLAB
clear

% Data files
years = 2015;
varnames = {'sst'};
% basedatapath = '~/Desktop/Geophysical_data/ERA-Interim/SST/'; baseoutpath = basedatapath; % FOR TESTS
basedatapath = '/Volumes/rap/martigalitapias/ERA-Interim/SST/'; % on taku-leifr
baseoutpath = '/Volumes/rap/martigalitapias/ERA-Interim/SST_46km/'; % on taku-leifr

% Grids
kmgrid2 = '46'; % 28, 37 or 46
gridpath = '~/Desktop/Grids_maps/grids/'; % TESTS and taku-leifr
grid1name = 'ERAi_45N.mat';
load([gridpath grid1name]); grid1 = [lon lat]; clear lon; clear lat
grid2name = ['grid' kmgrid2 'km_45N.txt'];
grid2 = dlmread([gridpath grid2name],' ');


%% Use index list to bin dmspt from grid1 to grid2

load(['indlist_ERAi_45N_to_' kmgrid2 'km.mat'])
tic
nrows = size(indlist,1);
ncols = length(varnames); % usual definition
% ncols = length(days); dataout = nan(nrows,ncols); % TESTS, alternative to compare consecutive images

toc
for vv = 1:length(varnames)
    for yyyy = years
        
        datapath = [basedatapath num2str(yyyy) '/']; % uncomment for taku-leifr
        flist = dir([datapath '/*.nc']); % TESTS
        
        for mmdd = 1:length(flist)
            
            disp(['Opening ' flist(mmdd).name])
            filepath = ([datapath flist(mmdd).name]);
            datain = ncread(filepath,varnames{vv});
            datain = datain(:) - 273.15;
            
            %             h = figure(); clf, set(gcf,'units','centimeters','position',[5 20 40 20])
            %             subplot(1,2,1)
            %             plot_sst_stereo(datain,8,grid1)
            
            % Binning/reprojection
            dataout = nan(nrows,1);
            for rr = 1:nrows
                dataout(rr,:) = nanmean(datain(indlist{rr,1}));
                if ~mod(rr,10000), rr, toc, end
            end
            tmp = [yyyy mmdd nanmean(dataout)];
            
            %             subplot(1,2,2)
            %             dataout(isnan(dataout)) = 10;
            %             plot_sst_stereo(dataout,8,grid2)
            
            % Write file
            dataout(isnan(dataout)) = -999;
            dummy1 = regexp(flist(mmdd).name,'\.nc','split');
            dummy2 = regexp(dummy1{1,1},'\ERAi_','split');
            datename = dummy2{1,2};
            %             outpath = sprintf('%s%s_%s.txt',baseoutpath,varnames{vv},datename); % TESTS
            outpath = sprintf('%s%i/%s_%s.txt',baseoutpath,yyyy,varnames{vv},datename); % taku-leifr
            dlmwrite(outpath,dataout,'\t');
            
            %             dlmwrite('./summary_MeanSstDay.txt',tmp,'-append');
            
            %             print(h,sprintf('~/Desktop/Heintzenberg_Leck_collaboration/SST_maps/SST_%s.png',datename),'-dpng');
            
        end
    end
end
toc
