% DEFINE THE PATH FOR FILES USED FOR GAP-FILLING

function fpaths_fill = find_paths_4fill(dirpath,sensor,sensorSST,years,ndperiod,ndays,period,kmgrid2,iy,ip)

% Minus 8D (year does not change for first years' first 8D period)
if ip == ndperiod(1) && iy <= years(1)
    fname_minus8D = sprintf('%c%cCLIM%03.0f_%0.0f%s.nc',sensor,sensorSST,ndperiod(end),ndays,period);
    fpath_minus8D = sprintf('%s%c%c_%0.0f%s_%skm/CLIM/%s',...
        dirpath,sensor,sensorSST,ndays,period,kmgrid2,fname_minus8D);
elseif ip == ndperiod(1)
    fname_minus8D = sprintf('%c%c%0.0f%03.0f_%0.0f%s.nc',sensor,sensorSST,iy-1,ndperiod(end),ndays,period);
    fpath_minus8D = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%s',...
        dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy-1,fname_minus8D);
else
    fname_minus8D = sprintf('%c%c%0.0f%03.0f_%0.0f%s.nc',sensor,sensorSST,iy,ip-8,ndays,period);
    fpath_minus8D = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%s',...
        dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy,fname_minus8D);
end
fpaths_fill.minus8D = fpath_minus8D;

% Plus 8D (year does not change for first years' first 8D period)
if ip == ndperiod(end) && iy >= years(end)
    fname_plus8D = sprintf('%c%cCLIM%03.0f_%0.0f%s.nc',sensor,sensorSST,ndperiod(1),ndays,period);
    fpath_plus8D = sprintf('%s%c%c_%0.0f%s_%skm/CLIM/%s',...
        dirpath,sensor,sensorSST,ndays,period,kmgrid2,fname_plus8D);
elseif ip == ndperiod(end)
    fname_plus8D = sprintf('%c%c%0.0f%03.0f_%0.0f%s.nc',sensor,sensorSST,iy+1,ndperiod(1),ndays,period);
    fpath_plus8D = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%s',...
        dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy+1,fname_plus8D);
else
    fname_plus8D = sprintf('%c%c%0.0f%03.0f_%0.0f%s.nc',sensor,sensorSST,iy,ip+8,ndays,period);
    fpath_plus8D = sprintf('%s%c%c_%0.0f%s_%skm/%0.0f/%s',...
        dirpath,sensor,sensorSST,ndays,period,kmgrid2,iy,fname_plus8D);
end
fpaths_fill.plus8D = fpath_plus8D;

% 8D CLIM
fname_8DCLIM = sprintf('%c%cCLIM%03.0f_%0.0f%s.nc',sensor,sensorSST,ip,ndays,period);
fpaths_fill.CLIM8D = sprintf('%s%c%c_%0.0f%s_%skm/CLIM/%s',...
    dirpath,sensor,sensorSST,ndays,period,kmgrid2,fname_8DCLIM);

% MOCLIM
eightD2month = {1:4; 5:7; 8:11; 12:14; 15:18; 19:22; 23:25; 26:29; 30:33; 34:37; 38:41; 42:46};
ii = 0;
while ii <= length(eightD2month)
    ii = ii + 1;
    tmp = eightD2month{ii};
    tmp = (tmp-1)*ndays+1;
    if sum(tmp == ip) == 1
        break
    end
end
fname_MOCLIM = sprintf('%c%cCLIM%02.0f_MONTH.nc',sensor,sensorSST,ii);
fpaths_fill.MOCLIM = sprintf('%s%c%c_MONTH_%skm/CLIM/%s',...
    dirpath,sensor,sensorSST,kmgrid2,fname_MOCLIM);


