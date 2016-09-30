% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Martí Galí Tàpias, 18 Feb 2016
% Improved 26 Sep 2016
tic

% May want to remove summary text files
% ! rm Feb2016_summary*

%% Some initial settings
varnameS  = {'chl_gsm' 'PP' 'dmspt_Asst_chloc' 'dmspt_Asst_chlgsm' 'dmspt_Asst_chlcota' 'Ice'}; % VERSION DMSPT
% varnameS  = {'dmspt_Asst_chlgsm'}; % VERSION STATS ONLY
years = 2003:2015; % normally 2003:2015
nyears = length(years);
period = 'MONTH';
ice_crit = 0.1;
kmgrid2 = '28'; % 28, 37 or 46 km macropixel size
outformat = 'netcdf'; % 'netcdf' r 'text'

%% Set file paths and name
dirpath = '/Volumes/rap/martigalitapias/binned_data/';
grid2path = '~/Desktop/Grids_maps/grids/grid';
outpath = '/Volumes/rap/martigalitapias/binned_data/';
sensor = 'A';
sensorSST = 'M';

%% Grid 2, conversion scheme
load(strcat('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_',kmgrid2,'km.mat'));
zbot2 = zbotmin;
grid2 = dlmread([grid2path kmgrid2 'km_45N.txt']);
lat2 = grid2(:,2);
npixels2 = length(lat2);

%% Binning


for ip = 1:12
    
    % Preallocate output: macropixel means (grid 2), n variables, ndays
    VARSOUT = nan(npixels2,length(varnameS));
    
    for iv = 1:length(varnameS)
        
        varname = varnameS{iv};
        
        % Preallocate data storage over nyears
        TMP2 = nan(npixels2,nyears);
        
        % Preallocate data storage over nyears statistics for grids 1 and 2
        npixels2MAR = nan(1,nyears);
        npixels2MAR65N = nan(1,nyears);
        
        for iy = years
            
            % Define first day of month for each year
            ndays = [31 28 31 30 31 30 31 31 30 31 30 31];
            if ~mod(iy,4)
                ndays(2) = ndays(2) + 1; % leap year
            end
            ndperiod = cumsum([1 ndays(1:(end-1))]);
            
            filename = sprintf('%c%c%0.0f%03.0f_%s.nc',sensor,sensorSST,iy,ndperiod(ip),period);
            file_test = ['grep ' filename ' list_MONTH.txt']; % file list in local folder
            status = system(file_test);
            if ~status % test file
                sprintf('File %s found',filename)
                filepath = sprintf('%s%c%c_%s_%skm/%0.0f/%s',...
                    dirpath,sensor,sensorSST,period,kmgrid2,iy,filename);
                ni=ncinfo(filepath);
                stra=char(ni.Variables.Name);
                var_test = strmatch(varname,stra,'exact');
                if ~isempty(var_test) % test variable
                    sprintf('Opening %s, variable %s',filename,varname)
                    var_grid2 = ncread(filepath,varname);
                    var_grid2(var_grid2==-999) = nan; % Note that ice data on grid 2 no longer has the "continent" or "non-covered area" flags
                    TMP2(:,id+1) = var_grid2;
                    
                    % Prepare marine pixel count for stats only for 1 DMSPt product
                    if strcmp('dmspt_Asst_chlgsm',varname)
                        % On grid2
                        Ice2 = ncread(filepath,'Ice');
                        npixels2MAR(id+1) = sum(zbot2<0 & Ice2<ice_crit); % npixels non-terrestrial with Ice<ice_crit
                        npixels2MAR65N(id+1) = sum(zbot2<0 & Ice2<ice_crit & lat2>=65); % npixels non-terrestrial with Ice<ice_crit and >65N
                    end
                end
            end
        end % loop on nyears
        
        % Average ndays. Note that -999 was converted to NaN before
        VARSOUT(:,iv) = nanmean(TMP2,2);
        VARSOUT(isnan(VARSOUT)) = -999;
        
        % Store complete stats only for 1 DMSPt product
        if strcmp('dmspt_Asst_chlgsm',varname)
            if ~status && ~isempty(var_test)
                % Summary statistics for all latitudes
                M2 = [ip summary_stats(TMP2,npixels2MAR)];
                % Repeat summary stats for latitudes >65
                TMP2(lat2<65,:) = [];
                M2_65 = [ip summary_stats(TMP2,npixels2MAR65N)];
            else
                M2 = [ip 0 0 0 0];
                M2_65 = [ip 0 0 0 0];
            end
            dlmwrite(sprintf('summary_%s_%sCLIM_%skm_%s.txt',varname,period,kmgrid2,date),M2,'-append')
            dlmwrite(sprintf('summary65N_%s_%sCLIM_%skm_%s.txt',varname,period,kmgrid2,date),M2_65,'-append')
        end
        
    end % loop on varnameS
    
    % Write netcdf or text file
    newvarnameS = varnameS;
    outname = sprintf('%s%c%c_%s_%skm/CLIM/%c%cCLIM%02.0f_%s.nc',...
        outpath,sensor,sensorSST,period,kmgrid2,sensor,sensorSST,ip,period);
    if ~isempty(VARSOUT)
        if strcmp(outformat,'netcdf')
            for iv = 1:length(newvarnameS)
                nccreate(outname,newvarnameS{iv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                ncwrite(outname,newvarnameS{iv},VARSOUT(:,iv));
            end
        elseif strcmp(outformat,'text')
            dlmwrite(outname,VARSOUT,'delimiter','\t','precision','%.4f');
        end
    end
    
end % loop on month periods

toc
