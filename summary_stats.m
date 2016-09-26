% SUMMARY STATS FOR (F)DMS(Pt) CHAIN DATA

function OUT = summary_stats(MAT_IN,npixelsMAR_in)

MAT_INmean = nanmean(MAT_IN,2);
cnt_data_ndays = mean(sum(~isnan(MAT_IN))); % mean coverage in n images with spatial binning gap filling
cnt_data_ndays_mean = sum(~isnan(MAT_INmean)); % mean coverage after temporal binning

OUT = [cnt_data_ndays...
        cnt_data_ndays/mean(npixelsMAR_in)... % mean fractional coverage prior to temporal binning
        cnt_data_ndays_mean...
        cnt_data_ndays_mean/mean(npixelsMAR_in)... % mean fractional coverage after temporal binning
        mean(sum(~isnan(MAT_IN),2))... % mean observations per pixel prior to binning
        nanmean(MAT_INmean)... % variable mean of means
        nanmean(MAT_IN(:))... % variable mean of all data
        nanmin(MAT_IN(:))... % variable min
        nanmax(MAT_IN(:))]; % variable max