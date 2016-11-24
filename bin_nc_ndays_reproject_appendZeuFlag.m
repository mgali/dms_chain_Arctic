% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Martí Galí Tàpias, 18 Feb 2016
% Improved 26 Sep 2016
tic

today = date;

%% Some initial settings
varnameS = {'Kd488'};
flagVarName = 'flag_dmspt_Asst';
varnameSout = {'Zeu' 'fstrat' 'fmix' 'fcocc'};
years = 2011:2015;
ndays = 8; % number of days averaged
ndperiod = 1 + ndays*(0:(365/ndays)); % defines first day of n-days period
ndperiod(ndperiod>366-floor(ndays/2)) = []; % remove periods starting too close to end of year relative to binning interval
period = 'D';
ice_crit = 0.1;
kmgrid2 = '28'; % 28, 37 or 46 km macropixel size

%% Set file paths and name
dirpath = '/Volumes/output-dev/Takuvik/Teledetection/Couleur/SORTIES/34_2_0/NOCLIM/'; % or 35_0_0 just for 2015
grid1path = '/Volumes/output-prod/Takuvik/Teledetection/All/Constant/';
grid2path = '~/Desktop/Grids_maps/grids/grid';
outpath = '/Volumes/rap/martigalitapias/binned_data/';
sensor = 'A';
sensorSST = 'M';

%% Grid 1
lat1 = ncread(strcat(grid1path,sensor,'45N.nc'),'lat');
npixels1 = length(lat1);

%% Grid 2, conversion scheme
load(strcat('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_',kmgrid2,'km.mat'));
zbot2 = zbotmin;
grid2 = dlmread([grid2path kmgrid2 'km_45N.txt']);
lat2 = grid2(:,2);
npixels2 = length(lat2);
load(['indlist_A45N_to_' kmgrid2 'km.mat']); % grid conversion scheme

%% Binning
for iy = years
    
    % Change version of PP files if year is 2015 or later
    if iy >= 2015
        dirpath = '/Volumes/output-dev/Takuvik/Teledetection/Couleur/SORTIES/35_0_0/NOCLIM/';
    end
    
    for ip = ndperiod
        
        % Preallocate output: macropixel means (grid 2), n variables, ndays
        VARSOUT = nan(npixels2,length(varnameSout));
        
        for iv = 1:length(varnameS)
            
            varname = varnameS{iv};
            
            % Preallocate data storage over ndays
            TMP2 = nan(npixels2,ndays);
            FLAGS_grid2 = nan(npixels2,3);
            TMP2S = nan(npixels2,ndays); % preallocate Strat fraction
            TMP2M = nan(npixels2,ndays); % preallocate Mix fraction
            TMP2C = nan(npixels2,ndays); % preallocate Cocc fraction
            
            nd = ndays-1;
            if length(ndperiod)>1 && ip==ndperiod(end)
                nd = 365-ip;
                if ~mod(iy,4)
                    nd = nd+1; % leap year
                end
            end
            for id = 0:nd
                filename = sprintf('%c%c%0.0f%03.0f_PP.nc',sensor,sensorSST,iy,ip+id);
                file_test = ['grep ' filename ' ' sensor sensorSST '_list.txt']; % file list in local folder
                status = system(file_test);
                if ~status
                    sprintf('File %s found',filename)
                    filepath = sprintf('%s%0.0f/%03.0f/%s',dirpath,iy,ip+id,filename);
                    ni=ncinfo(filepath);
                    stra=char(ni.Variables.Name);
                    var_test = strmatch(varname,stra,'exact');
                    if ~isempty(var_test) % test variable
                        
                        sprintf('Opening %s, variable %s',filename,varname)
                        
                        % Process Kd normally, calculate Zeu
                        var_grid1 = ncread(filepath,varname);
                        Ice1 = ncread(filepath,'Ice');
                        Ice1(Ice1==max(Ice1)) = 0; % assign non-covered marine pixels to Ice=0
                        Ice1(Ice1<0 | Ice1>1) = nan; % remove -999 values or values >1 (used as flag)var_grid1(var_grid1==-999) = nan;
                        var_grid1(Ice1 > ice_crit) = nan; % do not use ocean color data with Ice>0.1 (sub-pixel contamination not corrected)
                        var_grid1(var_grid1<=0) = nan; % remove negative chl values produced by gsm algos, remove -999 and potential zeros
                        var_grid1 = 4.6*((var_grid1).^-1); % calculate Zeu
                        var_grid2 = repbin_grid1_grid2(var_grid1,indlist); % reproject and bin from grid1 to grid2
                        TMP2(:,id+1) = var_grid2;
                        
                        % Process flags
                        flag_grid1 = ncread(filepath,flagVarName);
                        flag_grid1 = flag_grid1';
                        exclude = (Ice1 > ice_crit) | (isnan(var_grid1));
                        pre_flag = flag_grid1(~exclude,:);
                        pre_flag = cellstr(pre_flag); % cell array of strings
                        FLAGS_grid1 = nan(npixels1,3);
                        FLAGS_grid1(~exclude,1) = strcmp(pre_flag,'strat');
                        FLAGS_grid1(~exclude,2) = strcmp(pre_flag,'mix');
                        FLAGS_grid1(~exclude,3) = strcmp(pre_flag,'cocco');
                        
                        % Calculate bin Sum on grid2 of each algo flag
                        FLAGS_grid2(:,1) = repbinSum_grid1_grid2(FLAGS_grid1(:,1),indlist);
                        FLAGS_grid2(:,2) = repbinSum_grid1_grid2(FLAGS_grid1(:,2),indlist);
                        FLAGS_grid2(:,3) = repbinSum_grid1_grid2(FLAGS_grid1(:,3),indlist);
                        FLAGS_grid2 = FLAGS_grid2./(sum(FLAGS_grid2,2)*[1 1 1]); % calculate fraction of each algo used macropixel by macropixel.
                        TMP2S(:,id+1) = FLAGS_grid2(:,1);
                        TMP2M(:,id+1) = FLAGS_grid2(:,2);
                        TMP2C(:,id+1) = FLAGS_grid2(:,3);
                    end
                end
            end % loop on nd (number of days in each ndperiod)
        end % loop on varnameS
        
        % Average ndays. Note that -999 was converted to NaN before
        VARSOUT(:,1) = nanmean(TMP2,2);
        VARSOUT(:,2) = nanmean(TMP2S,2);
        VARSOUT(:,3) = nanmean(TMP2M,2);
        VARSOUT(:,4) = nanmean(TMP2C,2);
        VARSOUT(isnan(VARSOUT)) = -999;
        
        % Write netcdf file
        if ~isempty(VARSOUT)
            outname = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%c%c%0.0f%03.0f_%0.0f%s.nc',outpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
            for iv = 1:length(varnameSout)
                nccreate(outname,varnameSout{iv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                ncwrite(outname,varnameSout{iv},VARSOUT(:,iv));
            end
        end
    end % loop on nday periods
end % loop on years
toc

