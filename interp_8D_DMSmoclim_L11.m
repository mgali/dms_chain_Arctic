% Interpolate
clc, clear all
tic

load('~/Desktop/DMS_database_L10clim_BM04MLDclim/DMSmoclim_L11.mat')
MAT_GRID1 = DMS;
baseoutpath = '~/Desktop/Croft_collaboration/AM_8D_1deg/CLIM'; % on MBP
varnameout = 'dmsL11';

%% Add extremes (mid-December and mid-January climatology) to be able
% to actually interpolate afterwards, put average of Jan and Dec at boh temporal extremes
EXTREMES = nanmean(MAT_GRID1(:,:,[1 12]),3);
VAR1(:,:,1) = EXTREMES;
VAR1(:,:,14) = EXTREMES;
VAR1(:,:,2:13) = MAT_GRID1;

%% Create coordinates for meshgrid

x = 1:360;
y = 1:180;
z = [1 15:30:(30*12-15) 365];
[X,Y,Z] = meshgrid(x,y,z);
zq = 4:8:364;
[Xq,Yq,Zq] = meshgrid(x,y,zq);

%% Interpolate 3D and save both as *.mat and *.nc

VAR1q = interp3(X,Y,Z,VAR1,Xq,Yq,Zq);
DMS8D = VAR1q;
save(sprintf('%s/DMS8D.mat',baseoutpath),'DMS8D');

ndperiod = 1:8:361;

for ip = 1:length(ndperiod)
    
    sprintf('File for DOY = %i',ndperiod(ip))
    
    outpath = sprintf('%s/L11CLIM%03i_1deg_DMS.nc',baseoutpath,ndperiod(ip));
    nccreate(outpath,varnameout,'format','netcdf4','Dimensions',{'r' size(VAR1,1) 'c' size(VAR1,2)});
    ncwrite(outpath,varnameout,DMS8D(:,:,ip));
    
    % Plot to check
    h = figure(ip); imagesc(VAR1q(:,:,ip)), colorbar
    print(h,sprintf('%s/L11CLIM%03i_1deg_DMS.png',baseoutpath,ndperiod(ip)),'-dpng');
    close(h)
    
end

toc
