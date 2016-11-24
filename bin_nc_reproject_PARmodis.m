% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Martí Galí Tàpias, 18 Feb 2016
% Improved 26 Sep 2016
tic

today = date;

%% Some initial settings
varnameS = {'PAR'};
subvarnameS = {'par'};
years = 2003:2015; % 2003:2015 after test is done
ndays = 8; % number of days averaged
ndperiod = 1 + ndays*(0:(365/ndays)); % defines first day of n-days period
ndperiod(ndperiod>366-floor(ndays/2)) = []; % remove periods starting too close to end of year relative to binning interval
period = 'D';
kmgrid2 = '28'; % 28, 37 or 46 km macropixel size
outformat = 'netcdf'; % 'netcdf' or 'text'

%% Set file paths and name
dirpath = '/Volumes/taku-njall/MODISA/L3BIN/';
grid1path = '/Volumes/output-prod/Takuvik/Teledetection/All/Constant/';
grid2path = '~/Desktop/Grids_maps/grids/grid';
outpath = '/Volumes/rap/martigalitapias/binned_data/';
sensor = 'A';
sensorSST = 'M';

%% Grid 1
lat1 = ncread(strcat(grid1path,sensor,'45N.nc'),'lat');
zbot1 = ncread('/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/Bathymetre/Province_Zbot_MODISA_L3binV2.nc','Zbot');
npixels1 = length(lat1);

%% Grid 2, conversion scheme
load(strcat('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_',kmgrid2,'km.mat'));
zbot2 = zbotmin;
grid2 = dlmread([grid2path kmgrid2 'km_45N.txt']);
lat2 = grid2(:,2);
npixels2 = length(lat2);
load(['indlist_A45N_to_' kmgrid2 'km.mat']); % grid conversion scheme

%% Latitude on global grid, cutoff latitude "minlat"
lat0 = ncread(strcat(grid1path,'AWholeWorld.nc'),'lat');
minlat = 45;

%% Create file list
% system('ls /Volumes/taku-njall/MODISA/L3BIN/*/*/A*.L3b_8D_PAR.nc > ~/Desktop/Artic_DOSES/Binning/list_8D_PAR_MODISA.txt');


%% Binning
for iy = years
    
    for ip = ndperiod
        
        % Preallocate output: macropixel means (grid 2), n variables, ndays
        VARSOUT = nan(npixels2,length(varnameS));
        
        for iv = 1:length(varnameS)
            
            varname = varnameS{iv};
            subvarname = subvarnameS{iv};
            
            % Define last day of ndperiod
            nd = ndays-1;
            if length(ndperiod)>1 && ip==ndperiod(end)
                nd = 365-ip;
                if ~mod(iy,4)
                    nd = nd+1; % leap year
                end
            end
            
            filename = sprintf('%c%0.0f%03.0f%0.0f%03.0f.L3b_%0.0f%c_%s.nc',sensor,iy,ip,iy,ip+nd,ndays,period,varname);
            file_test = ['grep ' filename ' ' 'list_8D_PAR_MODISA.txt']; % file list in local folder
            status = system(file_test);
            if ~status
                sprintf('File %s found',filename)
                filepath = sprintf('%s%0.0f/%03.0f/%s',dirpath,iy,ip,filename);
                sprintf('Opening %s, variable %s',filename,subvarname)
                
                % Read NASA L3BIN file and subset >45N
                var = h5read(filepath, strcat('/level-3_binned_data/',subvarname));
                bl = h5read(filepath, '/level-3_binned_data/BinList');
                myvar = var.sum ./ double(bl.weights);
                var_grid1 = nan(length(lat0),1);
                bin_num = double(bl.bin_num);
                var_grid1(bin_num) = myvar;
                var_grid1(lat0 < minlat) = [];
                if length(var_grid1) == length(lat1)
                    var_grid2 = repbin_grid1_grid2(var_grid1,indlist); % reproject and bin from grid1 to grid2
                    VARSOUT(:,iv) = var_grid2; % Fill VARSOUT
                    VARSOUT(isnan(VARSOUT)) = -999; % Convert NaN back to -999
                else
                    error('Vector lengths do not match')
                end % test lengths match
            end % test file
        end % loop on varnameS
        
        % Write netcdf or text file
        if ~isempty(VARSOUT)
            if strcmp(outformat,'netcdf')
                outname = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%c%c%0.0f%03.0f_%0.0f%s.nc',outpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
                % outname =
                % sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%c%c%0.0f%03.0f_%0.0f%s_PARtest.nc',outpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period); 
                % % for TESTS
                for iv = 1:length(varnameS)
                    nccreate(outname,varnameS{iv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                    ncwrite(outname,varnameS{iv},VARSOUT(:,iv));
                end
            end
        end
        
    end % loop on nday periods
end % loop on years
toc

% bin_nc_ndaysCLIM % Uncomment if you want to run code to compute climatological means after the spatiotemporal binning

