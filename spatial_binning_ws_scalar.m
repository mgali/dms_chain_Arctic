% SPATIAL BINNING OF WIND SPEED DATA AND CALCULATION OF SCALAR WIND SPEED
clear

dateoutformat = 'YYYYDDD'; % YYYYDDD or YYYYMMDD
outformat = 'nc'; % 'nc' or 'txt'

% Grids
kmgrid2 = '28'; % 28, 37 or 46
gridpath = '~/Desktop/Grids_maps/grids/'; % TESTS and taku-leifr
grid1name = 'ERAi_45N.mat';
load([gridpath grid1name]); grid1 = [lon lat]; clear lon; clear lat
grid2name = ['grid' kmgrid2 'km_45N.txt'];
grid2 = dlmread([gridpath grid2name],' ');
npixels2 = size(grid2,1);

% Data files
years = 2007; % 1998:2015;
invarnames = {'u10' 'v10'};
outvarnames = {'ws'};
% basedatapath = '~/Desktop/Geophysical_data/ERA-Interim/'; baseoutpath = basedatapath; % FOR TESTS
basedatapath = '/Volumes/rap/martigalitapias/ERA-Interim/'; % on taku-leifr
baseoutpath = sprintf('/Volumes/rap/martigalitapias/ERA-Interim/WS_%skm/',kmgrid2); % on taku-leifr

%% Use index list to bin dmspt from grid1 to grid2

load(['indlist_ERAi_45N_to_' kmgrid2 'km.mat'])
tic
nrows = size(indlist,1);
ncols = length(outvarnames); % usual definition
% ncols = length(days); dataout = nan(nrows,ncols); % TESTS, alternative to compare consecutive images

toc
for vv = 1:length(outvarnames)
    for yyyy = years
        
        datapath1 = sprintf('%s%sm_DAY/%04d/',basedatapath,invarnames{1},yyyy);
        flist1 = dir(sprintf('%sERAi_%04d*.nc',datapath1,yyyy));
        
        datapath2 = sprintf('%s%sm_DAY/%04d/',basedatapath,invarnames{2},yyyy);
        flist2 = dir(sprintf('%sERAi_%04d*.nc',datapath2,yyyy));
        
        if length(flist1) ~= length(flist2)
            error('The number of u10 and v10 files does not match')
        else
            
            for mmdd = 1:length(flist1)
                
                disp(['Opening ' flist1(mmdd).name])
                filepath1 = ([datapath1 flist1(mmdd).name]);
                datain1 = ncread(filepath1,invarnames{1});
                
                disp(['Opening ' flist2(mmdd).name])
                filepath2 = ([datapath2 flist2(mmdd).name]);
                datain2 = ncread(filepath2,invarnames{2});
                
                % Calculate scalar wind speed here
                datain = sqrt(datain1.^2 + datain2.^2);
                datain = datain(:); % VECTOR FORMAT!
                
                % STOPEED CODING HERE. NEED TO CHECK PRIOR TO CONTINUE.
                % PLOT SIDE BY SIDE AND GO AHEAD!
                
                % Binning/reprojection
                dataout = repbin_grid1_grid2(datain,indlist);
                
                % Write file
                dataout(isnan(dataout)) = -999;
                dummy1 = regexp(flist1(mmdd).name,'\.nc','split');
                dummy2 = regexp(dummy1{1,1},'\ERAi_','split');
                datename = dummy2{1,2};
                
                if strcmp(dateoutformat,'YYYYMMDD')
                    outpath = sprintf('%s%s_DAY/%i/%s_%s.%s',baseoutpath,outvarnames{vv},yyyy,outvarnames{vv},datename,outformat); % taku-leifr
                elseif strcmp(dateoutformat,'YYYYDDD')
                    Y = str2double(datename(1:4));
                    M = str2double(datename(5:6));
                    D = str2double(datename(7:8));
                    DOY = yearday(D,M,Y,0,0,0);
                    outpath = sprintf('%s%s_DAY/%i/%s_%0.0f%03.0f.%s',baseoutpath,outvarnames{vv},Y,outvarnames{vv},Y,DOY,outformat); % taku-leifr
                end
                
                if strcmp(outformat,'txt')
                    dlmwrite(outpath,dataout,'\t');
                elseif strcmp(outformat,'nc')
                    nccreate(outpath,outvarnames{vv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                    ncwrite(outpath,outvarnames{vv},dataout);
                end
                
                % Store mean in separate file
                tmp = [yyyy mmdd nanmean(dataout)];
                dlmwrite('./summary_MeanWSDay.txt',tmp,'-append');
                
                % % Plot side by side datain and dataout, save
                % h = figure(); clf, set(gcf,'units','centimeters','position',[5 20 40 40])
                % subplot(2,2,1), plot_ws_stereo(double(datain1(:)),8,grid1)
                % subplot(2,2,2), plot_ws_stereo(double(datain2(:)),8,grid1)
                % subplot(2,2,3), plot_ws_stereo(datain,8,grid1)
                % subplot(2,2,4), plot_ws_stereo(dataout,8,grid2)
                % print(h,sprintf('%s%s_DAY/%i/%s_%s_ERAvs%skm.png',baseoutpath,outvarnames{vv},yyyy,outvarnames{vv},datename,kmgrid2),'-dpng');
                % close(h)
                
            end % days loop
        end % check number of files match
    end % years loop
end % variables loop (not really use here)
toc
