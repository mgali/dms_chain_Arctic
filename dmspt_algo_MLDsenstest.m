% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Improved 14 Oct 2016

tic

%% Some initial settings

varnameSout  = {'dmspt_chlgsm_mimoc' 'dmspt_chlgsm_godas'};
years = 2003:2015;
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

%% Set file paths and name

dirpath = '/Volumes/rap/martigalitapias/binned_data/';
sstpath = '/Volumes/rap/martigalitapias/ERA-Interim/SST_28km/';
godaspath = '/Volumes/rap/martigalitapias/GODAS/MLD_8D_28km/';
mimocpath = '/Volumes/rap/martigalitapias/MIMOC/MLD_8D_28km/CLIM/';
outpath = '/Volumes/rap/martigalitapias/binned_data/';
sensor = 'A';
sensorSST = 'M';

for iy = years
    for ip = ndperiod
        
        filename = sprintf('%c%c%0.0f%03.0f_%0.0f%s.nc',sensor,sensorSST,iy,ip,ndays,period);
        file_test = ['grep ' filename ' list_8D.txt']; % file list in local folder
        status = system(file_test);
        
        if ~status % test file exists
            
            sprintf('File %s found',filename)
            filepath = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%s',...
                dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy,filename);
            sprintf('Opening %s',filename)
            
            chl = ncread(filepath,'chl_gsm'); chl(chl==-999) = nan;
            zeu = ncread(filepath,'Zeu'); zeu(zeu==-999) = nan;
            sst = dlmread(sprintf('%s%0.0f/sst_%0.0f%03.0f.txt',sstpath,iy,iy,ip));
            
            for io = 1:length(varnameSout)
                
                if strcmp(varnameSout{io},'dmspt_chlgsm_mimoc')
                    mld = ncread(sprintf('%sMLD_%03.0f.nc',mimocpath,ip),'mld');
                elseif strcmp(varnameSout{io},'dmspt_chlgsm_godas')
                    mld = ncread(sprintf('%s%0.0f/MLD_%0.0f%03.0f.nc',godaspath,iy,iy,ip),'mld');
                end
                
                % DMSPt algo (call function)
                % No PIC data used due to less restrictive Chl calculation criterion in Takuvik chain
                var_out = dmspt_algorithm_RSE2015(chl,zeu,mld,sst,nan(size(chl)));
                
                % Save data
                outname = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%c%c%0.0f%03.0f_%0.0f%s_MLDsenstest.nc',...
                    outpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
                nccreate(outname,varnameSout{io},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                ncwrite(outname,varnameSout{io},var_out);
                
            end % loop on output variables
        end % test file
    end % loop on nday periods
end % loop on years

