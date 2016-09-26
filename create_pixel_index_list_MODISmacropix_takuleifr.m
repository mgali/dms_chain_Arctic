% SPATIAL BINNING TESTS WITH MATLAB
clear

addpath ~/Documents/PROGRAMES_RECURSOS_INFORMATICS/MATLAB/m_map/
macrogrids = {'28' '37' '46'};
zbot = ncread('/Volumes/output-prod/Takuvik/Teledetection/Couleur/SORTIES/Bathymetre/Province_Zbot_MODISA_L3binV2.nc','Zbot');

for ii = 1:length(macrogrids)
    
    % Grids
    
    kmgrid2 = macrogrids{ii}; % 28, 37 or 46
    
    gridpath = '~/Desktop/Grids_maps/grids/';
    grid1name = 'A45N.nc';
    grid1 = [ncread([gridpath grid1name],'lon') ncread([gridpath grid1name],'lat')];
    grid2name = ['grid' kmgrid2 'km_45N.txt'];
    grid2 = dlmread([gridpath grid2name]);
    load(strcat('~/Desktop/Geophysical_data/Bathymetry_gebco08/Grids_SW_MODIS_4gebco/gebco_08_',kmgrid2,'km.mat'));
    
    %% 1.1 Set new macropixel x and y length
    
    if strcmp(kmgrid2,'28'), factor = 6;
    elseif strcmp(kmgrid2,'37'), factor = 8;
    elseif strcmp(kmgrid2,'46'), factor = 10;
    end
    xlength = 4.64*factor;
    ylength = 4.64*factor;
    
    %% 1.2 Create index list linking grid1 to grid 2. Export array as *.mat
    % AND create index list linking grid2 to grid1 (REVERSE)
    
    allpixlist1 = 1:size(grid1,1);
    indlist_reverse = nan(size(allpixlist1));
    
    tic
    for i = 1:size(grid2,1)
        [lon2vm,lat2v] = find_polygon_vertexs(grid2(i,1),grid2(i,2),xlength,ylength);
        pix_match_ind = inpolygon(grid1(:,1),grid1(:,2),lon2vm,lat2v);
        indlist_reverse(pix_match_ind) = i;
        if zbotmin(i)<0
            indlist{i,1} = allpixlist1(pix_match_ind);
        else
            indlist{i,1} = nan;
        end
        if ~mod(i,100), i, kmgrid2, end
    end
    toc
    
    save(['indlist_A45N_to_' kmgrid2 'km.mat'],'indlist');
    save(['indlist_' kmgrid2 'km_to_A45N.mat'],'indlist_reverse');
    
    %% 1.3 Convert indlist to matrix. Export in txt format
    
    load(['indlist_A45N_to_' kmgrid2 'km.mat'])
    INDLIST = nan(size(indlist,1),factor^2);
    for i = 1:size(indlist,1)
        tmp = [];
        tmp = indlist{i,:};
        if length(tmp) >= factor^2
            tmp = tmp(1:factor^2);
            tmp(isnan(tmp)) = -999;
        elseif ~isempty(tmp)
            tmp(isnan(tmp)) = -999;
            tmp = [tmp -999*ones(1,factor^2-length(tmp))];
        end
        INDLIST(i,:) = tmp;
    end
    
    dlmwrite(['indlist_A45N_to_' kmgrid2 'km.txt'],INDLIST,' ');
    
    % Neet to remove since variable was not initialized
    clear indlist
    
end
