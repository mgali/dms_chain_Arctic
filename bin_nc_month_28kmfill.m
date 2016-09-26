% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Mart� Gal� T�pias, 18 Feb 2016
tic

% May want to remove summary text files
% ! rm Feb2016_summary*

varnameS  = {'chl_gsm' 'PP' 'dmspt_Asst_chloc' 'dmspt_Asst_chlgsm' 'dmspt_Asst_chlcota' 'Ice'}; % VERSION DMSPT
years = 2003:2015; % normally 2003:2015
ndays = [31 28 31 30 31 30 31 31 30 31 30 31];
period = 'MONTH';
ice_crit = 0.1;

%% Set file paths and name
dirpath = '/Volumes/output-dev/Takuvik/Teledetection/Couleur/SORTIES/34_2_0/NOCLIM/'; % or 35_0_0 just for 2015
gridpath = '/Volumes/output-prod/Takuvik/Teledetection/All/Constant/';
mldpath = '/Volumes/output-prod/Takuvik/Teledetection/All/Clim/MLD/';
outpath = '/Volumes/rap/martigalitapias/binned_data/';
sensor = 'A';
sensorSST = 'M';

%% Grid 1
lat1 = ncread(strcat(gridpath,sensor,'45N.nc'),'lat');
zbot1 = ncread('/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/Bathymetre/Province_Zbot_MODISA_L3binV2.nc','Zbot');
npixels1 = length(lat1);

%% Grid 2, conversion scheme
kmgrid2 = '28'; % 28, 37 or 46
load(strcat('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_',kmgrid2,'km.mat'));
zbot2 = zbotmin;
grid2 = dlmread(['~/Desktop/Grids_maps/grids/grid' kmgrid2 'km_45N.txt']);
lat2 = grid2(:,2);
npixels2 = length(lat2);
load(['indlist_A45N_to_' kmgrid2 'km.mat']); % grid conversion scheme
% load(['indlist_' kmgrid2 'km_to_A45N.mat']) % no reverse conversion (too
% long)

%% Binning
for iy = years
    
    % Define first day of month and number of days in month
    if ~mod(iy,4)
        ndays(2) = ndays(2) + 1; % leap year
    end
    % Change version of PP files if year is 2015 or later
    if iy >= 2015
        dirpath = '/Volumes/output-dev/Takuvik/Teledetection/Couleur/SORTIES/35_0_0/NOCLIM/';
    end
    ndperiod = cumsum([1 ndays(1:(end-1))]);
    mo = 0;
    for ip = ndperiod
        
        mo = mo + 1;
        
        % Preallocate output
        VARSOUT = nan(npixels2,length(varnameS));
        
        for iv = 1:length(varnameS)
            
            varname = varnameS{iv};
            
            % Preallocate data storage over ndays
            TMP1 = nan(npixels1,ndays(mo));
            TMP2 = nan(npixels2,ndays(mo));
            
            % Preallocate data storage over ndaysstatistics for grids 1 and 2
            npixels1MAR = nan(1,ndays(mo));
            npixels1MAR65N = nan(1,ndays(mo));
            npixels2MAR = nan(1,ndays(mo));
            npixels2MAR65N = nan(1,ndays(mo));
            
            for id = 0:(ndays(mo)-1)
                filename = sprintf('%c%c%0.0f%03.0f_PP.nc',sensor,sensorSST,iy,ip+id);
                file_test = ['grep ' filename ' ' sensor sensorSST '_list.txt']; % file list in local folder
                status = system(file_test);
                if ~status % test file
                    filepath = sprintf('%s%0.0f/%03.0f/%s',dirpath,iy,ip+id,filename);
                    ni=ncinfo(filepath);
                    stra=char(ni.Variables.Name);
                    var_test = strmatch(varname,stra,'exact');
                    if ~isempty(var_test) % test variable
                        var_grid1 = ncread(filepath,varname);
                        var_grid2 = repbin_grid1_grid2(var_grid1,indlist); % reproject and bin from grid1 to grid2
                        TMP1(:,id+1) = var_grid1;
                        TMP2(:,id+1) = var_grid2;
                        
                        % Prepare marine pixel count for stats only for 1 DMSPt product
                        if strcmp('dmspt_Asst_chlgsm',varname)
                            % On grid1
                            Ice1 = ncread(filepath,'Ice');
                            npixels1MAR(id+1) = sum(zbot1<0 & Ice1<ice_crit); % npixels non-terrestrial with Ice<ice_crit
                            npixels1MAR65N(id+1) = sum(zbot1<0 & Ice1<ice_crit & lat1>=65); % npixels non-terrestrial with Ice<ice_crit and >65N
                            % On grid2
                            Ice2 = repbin_grid1_grid2(Ice1,indlist);
                            npixels2MAR(id+1) = sum(zbot2<0 & Ice2<ice_crit); % npixels non-terrestrial with Ice<ice_crit
                            npixels2MAR65N(id+1) = sum(zbot2<0 & Ice2<ice_crit & lat2>=65); % npixels non-terrestrial with Ice<ice_crit and >65N
                        end
                    end
                end
            end % loop on nd (number of days in each ndperiod)
            
            
            % Store complete stats only for 1 DMSPt product
            if strcmp('dmspt_Asst_chlgsm',varname)
                if ~status && ~isempty(var_test)
                    % Summary statistics for all latitudes
                    M1 = [iy ip summary_stats(TMP1,npixels1MAR)];
                    M2 = [iy ip summary_stats(TMP2,npixels2MAR)];
                    % Repeat summary stats for latitudes >65
                    TMP1(lat1<65,:) = [];
                    TMP2(lat2<65,:) = [];
                    M1_65 = [iy ip summary_stats(TMP1,npixels1MAR65N)];
                    M2_65 = [iy ip summary_stats(TMP2,npixels2MAR65N)];
                else
                    M1 = [iy ip 0 0 0 0 0 nan nan nan nan];
                    M1_65 = [iy ip 0 0 0 0 0 nan nan nan nan];
                    M2 = [iy ip 0 0 0 0 0 nan nan nan nan];
                    M2_65 = [iy ip 0 0 0 0 0 nan nan nan nan];
                end
                dlmwrite(sprintf('summary_%s_%s_4km_%s.txt',varname,period,date),M1,'-append')
                dlmwrite(sprintf('summary_65N_%s_%s_4km_%s.txt',varname,period,date),M1_65,'-append')
                dlmwrite(sprintf('summary_%s_%s_28km_%s.txt',varname,period,date),M2,'-append')
                dlmwrite(sprintf('summary_65N_%s_%s_28km_%s.txt',varname,period,date),M2_65,'-append')
            end
            
            % Average ndays.
            if ~status && ~isempty(var_test)
                % TMP1(TMP1==-999) = nan; % no longer needed since -999 converted to nan before
                VARSOUT(:,iv) = nanmean(TMP2,2);
            end
            
        end % loop on varnameS
        
        % Write netcdf
        VARSOUT(isnan(VARSOUT)) = -999;
        newvarnameS = varnameS;
        ncoutname = sprintf('%s%c%c_%s_28km/%0.0f/%c%c%0.0f%03.0f_%s.nc',outpath,sensor,sensorSST,period,iy,sensor,sensorSST,iy,ip,period);
        if ~isempty(VARSOUT)
            for iv = 1:length(newvarnameS)
                nccreate(ncoutname,newvarnameS{iv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                ncwrite(ncoutname,newvarnameS{iv},VARSOUT(:,iv));
            end
        end
        
    end % loop on nday periods
end % loop on years
toc
