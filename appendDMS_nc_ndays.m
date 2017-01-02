% CALCULATE DMS FROM DMSPT AND PAR USING 2 MODELS WITH DIFFERENT PARAMETERS
% 9 DEc 2016

tic

%% Some initial settings

varnameS = {'dmspt_oc_filled' 'dmspt_gsm_filled' 'dmspt_cota_filled'};
splitnameS = regexp(varnameS,'_','split');
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

ice_crit = 0.1;

% DMS models
dmsmodels = {'dmsN' 'dmsB'};
dmseq = [-1.1842 0.5246 0.0206;
    -1.9240 1.0000 0.0196];

% Preallocate stats
COUNTS = nan(length(years)*length(ndperiod),length(varnameS)+2);
headcounts = {'y' '8dperiod' '8D_DMS_OC' '8D_DMS_GSM' '8D_DMS_COTA'};

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
                        
            for iv = 1:length(varnameS)
                
                varname = varnameS{iv};
                splitname = splitnameS{iv};
                var_in = ncread(filepath,varname); var_in(var_in<0) = nan;
                par = ncread(filepath,'PAR'); par(par==-999) = nan;
                Ice = ncread(filepath,'Ice'); Ice(Ice==-999) = nan;
                % Initial stats
                COUNTS(cc,[1 2 3+iv-1]) = [iy ip sum(~isnan(var_in) & Ice<=ice_crit)];
                
                % ///////// CALCULATE DMS ACCORDING TO TWO MODELS /////////
                for im = 1:length(dmsmodels)
                    
                    model = dmseq(im,:);
                    varnameout = sprintf('%s_%s_%s',dmsmodels{im},splitname{2},splitname{3});
                    var_out = 10.^(model(1) + model(2)*real(log10(var_in)) + model(3)*par);
                    var_out(Ice>ice_crit) = nan;
                    
                    % Save data
                    
                    var_out(isnan(var_out)) = -999;
                    
                    outname = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%c%c%0.0f%03.0f_%0.0f%s_DMS.nc',... % different file
                        outpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
                    nccreate(outname,varnameout,'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                    ncwrite(outname,varnameout,var_out);
                    
                end % loop on DMS models
            end % loop on variables
        end % test file
    end % loop on nday periods
    iy, toc
end % loop on years

save(strcat('stats_dms',date,'.mat'),'COUNTS','headcounts');
toc
