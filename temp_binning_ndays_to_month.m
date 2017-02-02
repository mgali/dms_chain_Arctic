% DETERMINE WEIGHTS OF N DAYS PERIODS TO CALCULATE MONTHLY CLIMATOLOGIES
% LEAP YEARS OMITTED

%% Calculate period weights for each month

% Define n days period and month
years = {'2004' '2005' '2006' '2007' '2008' '2009' '2010' '2011' '2012' '2013' '2014' '2015' 'CLIM'}; % '2003' 
% years = {'CLIM'};
ndays = 8; % number of days averaged
ndperiod = 1 + ndays*(0:(365/ndays)); % defines first day of n-days period
ndperiod(ndperiod>366-floor(ndays/2)) = []; % remove periods starting too close to end of year relative to binning interval

months = [1 31 28 31 30 31 30 31 31 30 31 30 31];
months = cumsum(months);

% Find weights of n days periods for each month
for mm = 1:(length(months)-1)
    mperiods = ndperiod(ndperiod>months(mm)-ndays & ndperiod<months(mm+1)); % [months(mm) months(mm+1)]
    pdays = ndays*ones(size(mperiods));
    pdays(1) = ndays + (mperiods(1) - months(mm));
    pdays(end) = months(mm+1) - mperiods(end); % pdays
    ww = pdays/sum(pdays);
    mp_nd_mo{mm} = mperiods;
    w_nd_mo{mm} = ww;
end

%% Define input variables and grids

varnameS = {'dmsN_oc_filled' 'dmsB_oc_filled' 'dmsN_gsm_filled' 'dmsB_gsm_filled' 'dmsN_cota_filled' 'dmsB_cota_filled'}; % filled dms
% varnameS = {'pic'}; % PIC, added to monthly clim through this express pathway (instead of computing means of monthly means from 4.6 km data)
% varnameS = {'PAR'}; % PAR, added to monthly clim through this express pathway (instead of computing means of monthly means from 4.6 km data)
period = 'D';
outperiod = 'MONTH';
kmgrid2 = '28'; % 28, 37 or 46 km macropixel size
outformat = 'netcdf'; % 'netcdf' r 'text'

%% Set file paths and name
% dirpath = '/Volumes/rap/martigalitapias/binned_data/'; % on taku-leifr
% outpath = '/Volumes/rap/martigalitapias/binned_data/'; % on taku-leifr
dirpath = '~/Desktop/Artic_DOSES/'; % on my MBP
outpath = '~/Desktop/Artic_DOSES/'; % on my MBP
grid2path = '~/Desktop/Grids_maps/grids/grid';
sensor = 'A';
sensorSST = 'M';
extraname = '_DMS'; % '_DMS', ''

%% Grid 2 (28 km is standard)
% load(strcat('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_',kmgrid2,'km.mat'));
grid2 = dlmread([grid2path kmgrid2 'km_45N.txt']);
npixels2 = size(grid2,1);

%% Binning

tic
for iy = 1:length(years)
    sprintf('8D to month binning for year %s',years{iy})
    for im = 1:12
        
        mp = mp_nd_mo{im};
        ww = w_nd_mo{im};
        VARSOUT = nan(npixels2,length(varnameS));
        
        for iv = 1:length(varnameS)
            
            varname = varnameS{iv};
            TMP = nan(npixels2,length(mp));
            
            for ip = 1:length(mp) % loop of n day periods overlapping with month
                filename = sprintf('%c%c%s%03i_%i%s%s.nc',sensor,sensorSST,years{iy},mp(ip),ndays,period,extraname);
                filepath = sprintf('%s%c%c_%i%s_%skm/%s/%s',dirpath,sensor,sensorSST,ndays,period,kmgrid2,years{iy},filename);
                TMP(:,ip) = ncread(filepath,varname);
            end
            TMP(TMP==-999) = nan;
            
            % Weighted average of n day periods. Note that -999 was converted to NaN before
            W = ones(npixels2,1)*ww;
            W(isnan(TMP)) = nan; % if no data, set weight to nan
            wrowsum = nansum(W,2);
            WROWSUM = wrowsum*ones(1,length(ww));
            WBIS = W./WROWSUM; % recalculate weights taking into account empty data
            WTMP = WBIS.*TMP;
            VARSOUT(:,iv) = nansum(WTMP,2); % note nansum of nans produces 0. Corrected in next line
            VARSOUT(isnan(VARSOUT) | VARSOUT==0) = -999;
            
        end % variables loop
        
        % Write netcdf or text file
        newvarnameS = varnameS;
        outname = sprintf('%s%c%c_%s_%skm/%s/%c%c%s%02i_%s%s.nc',...
            outpath,sensor,sensorSST,outperiod,kmgrid2,years{iy},sensor,sensorSST,years{iy},im,outperiod,extraname);
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
        
    end % months loop
end % years loop
toc
