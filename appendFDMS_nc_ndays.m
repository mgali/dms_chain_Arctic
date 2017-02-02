% CALCULATE DMS FROM DMSPT AND PAR USING 2 MODELS WITH DIFFERENT PARAMETERS
% 9 DEC 2016

tic

%% Some initial settings

varnameS = {'dmsB_oc_filled' 'dmsB_gsm_filled' 'dmsB_cota_filled'};
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

% Wind speed parameterization
param = 'W97';

%% Ice mask criterion and bottom depth mask

ice_crit = 0.1;
load(strcat('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_',kmgrid2,'km.mat'));
zbot2 = zbotmin;

%% Set file paths and name

dirpath = '/Volumes/rap/martigalitapias/binned_data';
meteopath = '/Volumes/rap/martigalitapias/ERA-Interim';
salwoapath = '~/Desktop/Geophysical_data/Salinity_WOA13';
outpath = '/Volumes/rap/martigalitapias/binned_data';
sensor = 'A';
sensorSST = 'M';

%% Process files

for iy = years
    for ip = ndperiod
        
        % Load ice, sst and sal_clim
        icepath = sprintf('%s/%c%c_%0.0f%s_%skm/%04i/%c%c%04i%03i_%i%s.nc',...
            dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
        Ice = ncread(icepath,'Ice'); Ice(Ice==-999) = nan;
        
        sstpath = sprintf('%s/SST_%skm/%04i/sst_%04i%03i.txt',...
            meteopath,kmgrid2,iy,iy,ip);
        sst = dlmread(sstpath); sst(sst==-999) = nan;
        
        salpath = sprintf('%s/SAL_%i%s_%skm/CLIM/SALCLIM%03i.nc',...
            salwoapath,ndays,period,kmgrid2,ip);
        sal = dlmread(salpath); sal(sal==-999) = nan;
        
        MAT_OUT = nan(npixels2,ndays);
        
        for iv = 1:length(varnameS)
            
            % Load nday dms
            varname = varnameS{iv};
            varpath = sprintf('%s/%c%c_%0.0f%s_%skm/%04i/%c%c%04i%03i_%i%s_DMS.nc',...
                dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
            var_in = ncread(varpath,varname); var_in(var_in<0) = nan;
            
            % Unobserved pixels: fill with mean of database measurements
            % where SZA at noon was >70 degrees, = 0.5 nM.
            var_in(zbot2<0 & Ice<=ice_crit & isnan(var_in)) = 0.5;
            
            % Define how many days searched within nday period
            pdays = ndays - 1;
            if ip == ndperiod(end)
                if mod(iy,4), pdays = 365 - ndperiod(end);
                else pdays = 366 - ndperiod(end);
                end
            end
            
            for id = 0:pdays
                
                % Load daily wind speed
                wspath = sprintf('%s/WS_%skm/%04i/ws_04i%03i.nc',...
                    meteopath,kmgrid2,iy,iy,ip);
                ws = dlmread(wspath); ws(ws==-999) = nan;
                
                varnameout = sprintf('%c%s','f',varname);
                var_out = fdms(var_in,ws,sst,sal,param);
                var_out(Ice>ice_crit) = nan;
                MAT_OUT(:,id+1) = var_out;
                
            end % loop on nday period days
            
            % Mean flux, taking into account that nan means zero
            out = nansum(MAT_OUT,2)/ndays;
            out(isnan(out)) = -999;
            
            % Save data
            outname = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%c%c%0.0f%03.0f_%0.0f%s_DMS.nc',... % different file
                outpath,sensor,sensorSST,ndays,period,kmgrid2,iy,sensor,sensorSST,iy,ip,ndays,period);
            nccreate(outname,varnameout,'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
            ncwrite(outname,varnameout,out);
            
        end % loop on variables
    end % loop on nday periods
    iy, toc
end % loop on years

toc
