% FILL 8D 28 KM DATA PROGRESSIVELY WITH MONTH AND MOCLIM DATA
% 24 Nov 2016

tic

%% Some initial settings

varnameS = {'dmspt_Asst_chloc' 'dmspt_Asst_chlgsm' 'dmspt_Asst_chlcota'};
varnameSout  = {'dmspt_oc_filled' 'dmspt_gsm_filled' 'dmspt_cota_filled'};
years = 2003:2015; % normally 2003:2015
ndays = 8; % number of days averaged
ndperiod = 1 + ndays*(0:(365/ndays)); % defines first day of n-days period
ndperiod(ndperiod>366-floor(ndays/2)) = []; % remove periods starting too close to end of year relative to binning interval
period = 'D';
kmgrid2 = '28'; % 28, 37 or 46 km macropixel size
if strcmp(kmgrid2,'28')
    npixels2 = 96668;
elseif strcmp(kmgrid2,'37')
    npixels2 = 54371;
elseif strcmp(kmgrid2,'46')
    npixels2 = 34799;
end

ice_crit_fill = 0.1;

% Preallocate stats
COUNTS = nan(length(years)*length(ndperiod),6);
headcounts = {'y' '8dperiod' '8D' 'plus24D' 'plus8DCLIM' 'plusMOCLIM'};
var4counts = 'dmspt_Asst_chlgsm';

%% Set file paths and name

dirpath = '/Volumes/rap/martigalitapias/binned_data/';
outpath = '/Volumes/rap/martigalitapias/binned_data/';
sensor = 'A';
sensorSST = 'M';


%% Process files
cc = 0;
for iy = years
    for ip = ndperiod
        
        filename = sprintf('%c%c%0.0f%03.0f_%0.0f%s.nc',sensor,sensorSST,iy,ip,ndays,period);
        file_test = ['grep ' filename ' list_8D.txt']; % file list in local folder
        status = system(file_test);
        cc = cc + 1;
        
        if ~status % test file exists
            
            sprintf('File %s found',filename)
            filepath = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%s',...
                dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy,filename);
            sprintf('Opening %s',filename)
            
            % Define file paths of ±8D, 8DCLIM and monthly clim data
            fpaths_fill = find_paths_4fill(dirpath,sensor,sensorSST,years,ndperiod,ndays,period,kmgrid2,iy,ip);
            
            for iv = 1:length(varnameS)
                
                varname = varnameS{iv};
                var_out = ncread(filepath,varname); var_out(var_out==-999) = nan;
                Ice = ncread(filepath,'Ice'); Ice(Ice==-999) = nan;
                % Initial stats
                if strcmp(varname,var4counts), COUNTS(cc,1:3) = [iy ip sum(~isnan(var_out))]; end
                
                % FILL WITH TEMPORALLY BINNED DATA
                var_fillminus8D = ncread(fpaths_fill.minus8D,varname); var_fillminus8D(var_fillminus8D==-999) = nan;
                var_fillplus8D = ncread(fpaths_fill.plus8D,varname); var_fillplus8D(var_fillplus8D==-999) = nan;
                var_fillpm8Dt = ([var_fillminus8D var_fillplus8D])';
                var_fill24D = (nanmean(var_fillpm8Dt))';
                var_fill8DC = ncread(fpaths_fill.CLIM8D,varname); var_fill8DC(var_fill8DC==-999) = nan;
                var_fillMC = ncread(fpaths_fill.MOCLIM,varname); var_fillMC(var_fillMC==-999) = nan;
                
                % First fill step: prior and posterior 8D
                crit_fill24D = isnan(var_out) & ~isnan(var_fill24D) & Ice < ice_crit_fill;
                var_out(crit_fill24D == 1) = var_fill24D(crit_fill24D == 1);
                % First fill step STATS
                if strcmp(varname,var4counts), COUNTS(cc,4) = sum(~isnan(var_out)); end
                
                % Second fill step: climatological 8D
                crit_fill8DC = isnan(var_out) & ~isnan(var_fill8DC) & Ice < ice_crit_fill;
                var_out(crit_fill8DC == 1) = var_fill8DC(crit_fill8DC == 1);
                % Second fill step STATS
                if strcmp(varname,var4counts), COUNTS(cc,5) = sum(~isnan(var_out)); end
                
                % Third fill step: climatological month
                crit_fillMC = isnan(var_out) & ~isnan(var_fillMC) & Ice < ice_crit_fill;
                var_out(crit_fillMC == 1) = var_fillMC(crit_fillMC == 1);
                % Third fill step STATS
                if strcmp(varname,var4counts), COUNTS(cc,6) = sum(~isnan(var_out)); end
                
                
                % Save data
                
                var_out(isnan(var_out)) = -999;
                
                outname = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%c%c%0.0f%03.0f_%0.0f%s.nc',... % for tests
                    outpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
                nccreate(outname,varnameSout{iv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                ncwrite(outname,varnameSout{iv},var_out);
                
            end % loop on variables
        end % test file
    end % loop on nday periods
end % loop on years

save(strcat('stats_fill_dmspt',date,'.mat'),'COUNTS','headcounts','var4counts');

