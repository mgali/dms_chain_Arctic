% BIN DAILY IMAGES TO SPECIFIED PERIOD (NDAYS)
% Martí Galí Tàpias, 18 Feb 2016
% Improved 26 Sep 2016, 14 Apr 2017
tic
today = date;

%% Initial settings TO EDIT: variable names and corresponding paths

years = 1998:2007; % normally 2003:2016 for MODIS and 1998:2007 for SeaWiFS
ndays = 8; % number of days averaged
ndperiod = 1 + ndays*(0:(365/ndays)); % vector with first day of n-days periods
ndperiod(ndperiod>366-floor(ndays/2)) = []; % remove periods starting too close to end of year relative to binning interval
period = 'MONTH';
kmgrid2 = '28'; % 28, 37 or 46 km macropixel size
outformat = 'netcdf'; % 'netcdf' or 'text'
sensor = 'S'; % A for MODIS-Aqua, S for SeaWiFS
sensorATM = 'I'; % M for MODIS atmosphere, I for ISCCP
ice_crit = 0.1; % ice concentration above which data are not used

% ---------------- Variable groups and corresponding paths ----------------
filetype = 'dmspt'; % 'PP', 'dmspt', 'pic'

if strcmp(filetype,'PP')
    varnameS = {'chl_gsm' 'chl_oc' 'PP' 'Ice' 'aCDOM_412'};
    var4stats = 'chl_gsm';
    dirpath = '/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/35_3_2/NOCLIM';
    icedirpath = dirpath;
else
    icedirpath = '/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/35_3_2/NOCLIM';
    if strcmp(filetype,'dmspt')
        varnameS = {'dmspt_Asst_chlgsm'};% 'dmspt_Asst_chloc' 'dmspt_Asst_chlcota'};
        var4stats = 'dmspt_Asst_chlgsm';
        dirpath = '/Volumes/scratch/martigalitapias/35_3_2/NOCLIM';
    elseif strcmp(filetype,'pic')
        varnameS = {'pic'};
        dirpath = '/Volumes/output-prod/Takuvik/Teledetection/Products/pic/0_0_1_0';
        sensorATM = '';
    end
end

% File lists: Remove pre-existing lists and create up-to-date ones
% Lists created in loop (append), otherwise list is too long for bash.
% Comment if existing lists are OK.
dummy01 = system(sprintf('rm %c%c_list_PP_%s.txt',sensor,sensorATM,today));
if ~strcmp(filetype,'PP')
    dummy02 = system(sprintf('rm %c%c_list_%s_%s.txt',sensor,sensorATM,filetype,today));
end
for iy = years
    dummy1 = system(sprintf('ls %s/%04i/*/%c%c*_PP.nc >> %c%c_list_PP_%s.txt',icedirpath,iy,sensor,sensorATM,sensor,sensorATM,today));
    if ~strcmp(filetype,'PP')
        dummy2 = system(sprintf('ls %s/%04i/*/%c%c*_%s.nc >> %c%c_list_%s_%s.txt',dirpath,iy,sensor,sensorATM,filetype,sensor,sensorATM,filetype,today));
    end
end
toc, sprintf('New file lists created')

% ------------------------------ Common paths -----------------------------
grid1path = '/Volumes/output-prod/Takuvik/Teledetection/Grid/trunk/201510151636';
% grid1path = '~/Desktop/Grids_maps/grids/'; % Alternative (Desktop of my MBP or taku-leifr)
grid2path = '~/Desktop/Grids_maps/grids'; % Currently on my Desktop (on taku-leifr)
user = 'martigalitapias';
outpath = sprintf('/Volumes/rap/%s/binned_data',user);

% ---------------- Load grid 1 and corresponding bathymetry ---------------
lat1 = ncread(sprintf('%s/%c45N.nc',grid1path,sensor),'lat');
if sensor=='A'
    zbot1 = ncread('/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/Bathymetre/Province_Zbot_MODISA_L3binV2.nc','Zbot');
elseif sensor=='S'
    zbot1 = ncread('/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/Bathymetre/Province_Zbot_SeaWiFS_L3binV2.nc','Zbot');
end
npixels1 = length(lat1);

% ---------------- Load grid 2 and corresponding bathymetry ---------------
load(sprintf('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_%skm.mat',kmgrid2));
zbot2 = zbotmin; % may want tho choose zbotmean, more restrictive pixel inclusion
grid2 = dlmread(sprintf('%s/grid%skm_45N.txt',grid2path,kmgrid2));
lat2 = grid2(:,2);
npixels2 = length(lat2);

% --------------- Load grid1 -> grid2 conversion scheme -------------------
load(sprintf('indlist_%c45N_to_%skm.mat',sensor,kmgrid2)); % in current dir


%% Binning
for iy = years
    
    ndays = [31 28 31 30 31 30 31 31 30 31 30 31];
    % Define first day of month and number of days in month
    if ~mod(iy,4)
        ndays(2) = ndays(2) + 1; % leap year
    end
    
    ndperiod = cumsum([1 ndays(1:(end-1))]);
    mo = 0;
    for ip = ndperiod
        
        mo = mo + 1;
        
        % Preallocate output: macropixel means (grid 2), n variables, ndays
        VARSOUT = nan(npixels2,length(varnameS));
        
        for iv = 1:length(varnameS)
            
            varname = varnameS{iv};
            
            % Preallocate data storage over ndays
            TMP1 = nan(npixels1,ndays(mo));
            TMP2 = nan(npixels2,ndays(mo));
            
            % Preallocate data storage over ndays statistics for grids 1 and 2
            npixels1MAR = nan(1,ndays(mo));
            npixels1MAR65N = nan(1,ndays(mo));
            npixels2MAR = nan(1,ndays(mo));
            npixels2MAR65N = nan(1,ndays(mo));
            
            for id = 0:(ndays(mo)-1)
                filename = sprintf('%c%c%04i%03i_%s.nc',sensor,sensorATM,iy,ip+id,filetype);
                file_test = sprintf('grep %s %c%c_list_%s_%s.txt',filename,sensor,sensorATM,filetype,today); % test against file list in local folder
                status = system(file_test);
                icefilename = filename;
                status_ice = status;
                if ~strcmp(filetype,'PP')
                    icefilename = sprintf('%c%c%04i%03i_PP.nc',sensor,sensorATM,iy,ip+id);
                    ice_file_test = sprintf('grep %s %c%c_list_PP_%s.txt',icefilename,sensor,sensorATM,today); % test against file list in local folder
                    status_ice = system(ice_file_test);
                end
                if ~status && ~status_ice
                    sprintf('Files found')
                    filepath = sprintf('%s/%04i/%0i/%s',dirpath,iy,ip+id,filename);
                    icepath = sprintf('%s/%04i/%0i/%s',icedirpath,iy,ip+id,icefilename);
                    ni=ncinfo(filepath);
                    stra=char(ni.Variables.Name);
                    var_test = strmatch(varname,stra,'exact');
                    if ~isempty(var_test) % test variable
                        sprintf('Opening %s, variable %s',filename,varname)
                        var_grid1 = ncread(filepath,varname);
                        var_grid1(var_grid1==-999) = nan;
                        if strcmp('Ice',varname)
                            var_grid1(var_grid1==max(var_grid1)) = 0;
                            var_grid1(var_grid1<0 | var_grid1>1) = nan;
                        else
                            Ice1 = ncread(icepath,'Ice');
                            Ice1(Ice1==max(Ice1)) = 0; % assign non-covered marine pixels to Ice=0
                            Ice1(Ice1<0 | Ice1>1) = nan; % remove -999 values or values >1 (used as flag)
                            var_grid1(Ice1 > ice_crit) = nan; % do not use ocean color data with Ice>0.1 (sub-pixel contamination not corrected)
                            if ~(strcmp('Nsst',varname) || strcmp('sst_avhrr_4ql',varname))
                                var_grid1(var_grid1<0) = nan; % remove negative chl values produced by gsm algos
                            end
                        end
                        var_grid2 = repbin_grid1_grid2(var_grid1,indlist); % reproject and bin from grid1 to grid2
                        TMP1(:,id+1) = var_grid1;
                        TMP2(:,id+1) = var_grid2;
                        
                        % Prepare marine pixel count for stats only for 1 DMSPt product
                        if strcmp(var4stats,varname)
                            % On grid1
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
            
            % Average ndays. Note that -999 was converted to NaN before
            VARSOUT(:,iv) = nanmean(TMP2,2);
            VARSOUT(isnan(VARSOUT)) = -999;
            
            % Store complete stats only for chosen variable (1 is enough!)
            if strcmp(var4stats,varname)
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
                dlmwrite(sprintf('summary_%s_%s_%cgrid_%s.txt',varname,period,sensor,today),M1,'-append')
                dlmwrite(sprintf('summary65N_%s_%s_%cgrid_%s.txt',varname,period,sensor,today),M1_65,'-append')
                dlmwrite(sprintf('summary_%s_%s_%skm_%s.txt',varname,period,kmgrid2,today),M2,'-append')
                dlmwrite(sprintf('summary65N_%s_%s_%skm_%s.txt',varname,period,kmgrid2,today),M2_65,'-append')
            end
            
        end % loop on varnameS
        
        % Write netcdf or text file
        if ~isempty(VARSOUT)
            if strcmp(outformat,'netcdf')
                outname = sprintf('%s/%c%c_%s_%skm/%0.0f/%c%c%0.0f%03.0f_%s.nc',outpath,sensor,sensorATM,period,kmgrid2,iy,sensor,sensorATM,iy,ip,period);
                for iv = 1:length(varnameS)
                    nccreate(outname,varnameS{iv},'format','netcdf4','Dimensions',{'r' npixels2 'c' 1});
                    ncwrite(outname,varnameS{iv},VARSOUT(:,iv));
                end
            elseif strcmp(outformat,'text')
                outname = sprintf('%s/%c%c_%s_%skm/%0.0f/%c%c%0.0f%03.0f_%s.txt',outpath,sensor,sensorATM,period,kmgrid2,iy,sensor,sensorATM,iy,ip,period);
                dlmwrite(outname,VARSOUT,'delimiter','\t','precision','%.4f');
            end
        end
        
    end % loop on nday periods
end % loop on years
toc

