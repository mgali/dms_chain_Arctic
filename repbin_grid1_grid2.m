% SPATIAL BINNING FROM GRID1 (HIGHER RESOLUTION) TO GRID2

function [dataout] = repbin_grid1_grid2(datain,indexlist_array)

nrows = length(indexlist_array);
dataout = nan(length(indexlist_array),1);
for rr = 1:nrows
    match = indexlist_array{rr,1};
    if isnan(match(1))
        continue
    else
    dataout(rr) = nanmean(datain(match));
    end
end
