% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Improved 14 Oct 2016

function extract_regionStats_nc_ndays(num_region)

tic

%% Some initial settings

regions = {'Norway' 'Greenland' 'IcelandS' 'GulfAlaska' 'BeringS' 'Chukchi' 'BarentsS' 'Labrador' 'BaffinBay'};
% varnameS  = {'dmspt_Asst_chlgsm' 'dmspt_Asst_chloc' 'dmspt_Asst_chlcota' 'Zeu' 'fstrat' 'fmix' 'fcocc' 'chl_gsm' 'PP'}; % for main chain
varnameS  = {'dmspt_chlgsm_mimoc' 'dmspt_chlgsm_godas'};
years = 2003:2015;
ndays = 8; % number of days averaged
ndperiod = 1 + ndays*(0:(365/ndays)); % defines first day of n-days period
ndperiod(ndperiod>366-floor(ndays/2)) = []; % remove periods starting too close to end of year relative to binning interval
period = 'D';
kmgrid2 = '28'; % 28, 37 or 46 km macropixel size

%% Set file paths and name
dirpath = '/Volumes/rap/martigalitapias/binned_data/';
outpath = '~/Desktop/Artic_DOSES/Arctic_DOSES_timeseries/';
sensor = 'A';
sensorSST = 'M';

for ir = num_region
    % //////////////////////////////////////////////////////
    region = regions{ir};
    % //////////////////////////////////////////////////////
    
    %% Load region indexs
    load(strcat('~/Desktop/Artic_DOSES/Binning/indlist_',region,'_',kmgrid2,'km.mat'));
    
    %% Extraction of statistics
    
    statslist = {'mean' 'std' 'median' 'Q25' 'Q75' 'Q05' 'Q95' 'count'};
    
    % Preallocate output: macropixel means (grid 2), n variables, ndays
    STATSOUT = nan(length(ndperiod)*length(years),length(statslist),length(varnameS));
    
    rr = 0;
    
    for iy = years
        for ip = ndperiod
            
            rr = rr + 1;
            
            for iv = 1:length(varnameS)
                varname = varnameS{iv};
                % filename = sprintf('%c%c%0.0f%03.0f_%0.0f%s.nc',sensor,sensorSST,iy,ip,ndays,period);
                filename = sprintf('%c%c%0.0f%03.0f_%0.0f%s_MLDsenstest.nc',sensor,sensorSST,iy,ip,ndays,period);
%                 file_test = ['grep ' filename ' list_8D.txt']; % file list in local folder
%                 status = system(file_test);
%                 if ~status
                    sprintf('File %s found',filename)
                    filepath = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%s',...
                        dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy,filename);
                    ni=ncinfo(filepath);
                    stra=char(ni.Variables.Name);
                    var_test = strmatch(varname,stra,'exact');
                    if ~isempty(var_test) % test variable
                        sprintf('Opening %s, variable %s',filename,varname)
                        var_grid2 = ncread(filepath,varname);
                        regvar = var_grid2(regindlist);
                        regvar(regvar==-999) = [];
                        if ~isempty(regvar)
                            % Fill stats matrix
                            STATSOUT(rr,:,iv) = [nanmean(regvar)...
                                nanstd(regvar)...
                                nanmedian(regvar)...
                                quantile(regvar,[.25 .75 .05 .95])...
                                length(regvar)];
                        end
%                     end
                end
            end % loop on varnameS
        end % loop on nday periods
    end % loop on years
    
    toc
    
    stats_years = years;
    stats_ndperiod = ndperiod;
    stats_varname = varnameS;
    % save([outpath 'regStats_' region '_' num2str(ndays) period '_' kmgrid2 'km.mat'],'stats_years','stats_ndperiod','STATSOUT','statslist','stats_varname')
    save([outpath 'regStats_MLDsenstest_' region '_' num2str(ndays) period '_' kmgrid2 'km.mat'],'stats_years','stats_ndperiod','STATSOUT','statslist','stats_varname')
    
end

