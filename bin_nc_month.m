% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Martí Galí Tàpias, 18 Feb 2016
tic

% May want to remove some summary files
% ! rm Feb2016_summary*

% varnameS = {'dmspt_Asst_chlgsm'}; % TEST
% varnameS = {'dmspt_Asst_chlgsm' 'CF_mean'}; % VERSION FOR JOST HEINTZENBERG
% varnameS = {'Ice' 'chl_gsm' 'chl_gsm_mustapha' 'chl_cota' 'PP' 'PAR_cloud' 'CF_mean'}; % VERSION FOR SOPHIE
varnameS  = {'chl_gsm' 'PP' 'dmspt_Asst_chloc' 'dmspt_Asst_chlgsm' 'dmspt_Asst_chlcota' 'Ice' 'PAR_cloud'}; % VERSION DMSPT
version = ''; % Normally '', may want to add version number to end of file for tests (eg _34_2_0)
years = 2015:2015;
% years = 2003:2015;
ndays = [31 28 31 30 31 30 31 31 30 31 30 31];
period = 'MONTH';

% Set file paths and name
dirpath = '/Volumes/output-dev/Takuvik/Teledetection/Couleur/SORTIES/35_0_0/NOCLIM/'; % or 34_2_0
gridpath = '/Volumes/output-prod/Takuvik/Teledetection/All/Constant/';
mldpath = '/Volumes/output-prod/Takuvik/Teledetection/All/Clim/MLD/';
outpath = '/Volumes/rap/martigalitapias/binned_data/';
sensor = 'A';
sensorSST = 'M';
lat = ncread(strcat(gridpath,sensor,'45N.nc'),'lat');
lon = ncread(strcat(gridpath,sensor,'45N.nc'),'lon');
zbot = ncread('/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/Bathymetre/Province_Zbot_MODISA_L3binV2.nc','Zbot');
npixels = length(lat);
plist = 1:npixels;

for iy = years
    
    % Define first day of month and number of days in month
    if ~mod(iy,4)
        ndays(2) = ndays(2) + 1; % leap year
    end
    ndperiod = cumsum([1 ndays(1:(end-1))]);
    mo = 0;
    for ip = ndperiod
        
        mo = mo + 1;
        VARSOUT = nan(npixels,length(varnameS)); % preallocate
        
        for iv = 1:length(varnameS)
            
            varname = varnameS{iv};
            TMP = nan(npixels,ndays(mo)); % preallocate
            npixelsMAR = nan(1,ndays(mo));
            npixelsMAR65N = nan(1,ndays(mo));
            
            for id = 0:(ndays(mo)-1)
                filename = sprintf('%c%c%0.0f%03.0f_PP%s.nc',sensor,sensorSST,iy,ip+id,version);
                file_test = ['grep ' filename ' ' sensor sensorSST '_list.txt']; % file list in local folder
                status = system(file_test);
                if ~status
                    filepath = sprintf('%s%0.0f/%03.0f/%s',dirpath,iy,ip+id,filename);
                    ni=ncinfo(filepath);
                    stra=char(ni.Variables.Name);
                    var_test = strmatch(varname,stra,'exact');
                    if ~isempty(var_test)
                        TMP(:,id+1) = ncread(filepath,varname);
                        Ice = ncread(filepath,'Ice');
                        npixelsMAR(id+1) = sum(zbot<=0 & Ice<0.15); % npixels non-terrestrial with Ice<0.15
                        npixelsMAR65N(id+1) = sum(zbot<=0 & Ice<0.15 & lat>=65); % npixels non-terrestrial with Ice<0.15 and >65N                        
                        % disp(['Opening ' filename]); toc
                    end
                end
            end % loop on nd (number of days in each ndperiod)
            
            if ~status && ~isempty(var_test)
                % Average ndays. Summary statistics
                TMP(TMP==-999) = nan;
                TMPmean = nanmean(TMP,2);
                VARSOUT(:,iv) = TMPmean;
                cnt_data_ndays = mean(sum(~isnan(TMP))); % mean coverage in n images
                frac_data_ndays = cnt_data_ndays/mean(npixelsMAR);
                cnt_data_ndays_mean = sum(~isnan(TMPmean)); % mean coverage after temporal binning
                frac_data_ndays_mean = cnt_data_ndays_mean/mean(npixelsMAR);
                frac_2ormore = sum(sum(~isnan(TMP),2)>=2)/mean(npixelsMAR);
                mn = nanmean(TMPmean);
                md = nanmean(TMP(:));
                mi = nanmin(TMP(:));
                ma = nanmax(TMP(:));
                M = [iy ip cnt_data_ndays frac_data_ndays cnt_data_ndays_mean frac_data_ndays_mean frac_2ormore mn md mi ma];
                
                % Repeat summary stats for latitudes >65
                TMP(lat<65,:) = [];
                TMPmean = nanmean(TMP,2);
                cnt_data_ndays = mean(sum(~isnan(TMP))); % mean coverage in n images
                frac_data_ndays = cnt_data_ndays/mean(npixelsMAR65N);
                cnt_data_ndays_mean = sum(~isnan(TMPmean)); % mean coverage after temporal binning
                frac_data_ndays_mean = cnt_data_ndays_mean/mean(npixelsMAR65N);
                frac_2ormore = sum(sum(~isnan(TMP),2)>=2)/mean(npixelsMAR65N);
                mn = nanmean(TMPmean);
                md = nanmean(TMP(:));
                mi = nanmin(TMP(:));
                ma = nanmax(TMP(:));
                M65 = [iy ip cnt_data_ndays frac_data_ndays cnt_data_ndays_mean frac_data_ndays_mean frac_2ormore mn md mi ma];
            else
                M = [iy ip 0 0 0 0 0 nan nan nan nan];
                M65 = [iy ip 0 0 0 0 0 nan nan nan nan];
            end
            dlmwrite(sprintf('summary_%s_%s_4km.txt',varname,period),M,'-append')
            dlmwrite(sprintf('summary_65N_%s_%s_4km.txt',varname,period),M65,'-append')
            
        end % loop on varnameS
        
        % % Uncomment if lat lon added and empty rows not saved
        % VARSOUT = [lat lon VARSOUT];
        % VARSOUT(zbot>=0,:) = [];
        % VARSOUT(sum(isnan(VARSOUT),2)==length(varnameS),:) = [];
        % newvarnameS = ['lat' 'lon' varnameS];
        
        % Write netcdf
        VARSOUT(isnan(VARSOUT)) = -999;
        newvarnameS = varnameS;
        ncoutname = sprintf('%s%c%c_%s_4km/%0.0f/%c%c%0.0f%03.0f_%s.nc',outpath,sensor,sensorSST,period,iy,sensor,sensorSST,iy,ip,period);
        if ~isempty(VARSOUT)
            for iv = 1:length(newvarnameS)
                nccreate(ncoutname,newvarnameS{iv},'format','netcdf4','Dimensions',{'r' npixels 'c' 1});
                ncwrite(ncoutname,newvarnameS{iv},VARSOUT(:,iv));
            end
        end
        
    end % loop on nday periods
end % loop on years
toc
