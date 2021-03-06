% SUMMARY STATS FOR (F)DMS(Pt) CHAIN DATA

function OUT = summary_statsCLIM(MAT_IN,npixelsMAR_in)

MAT_INmean = nanmean(MAT_IN,2); % a npixels x 1 vector

cnt_data_ndays = sum(~isnan(MAT_IN)); % a 1 x ndays vector
cnt_data_ndays_mean = sum(~isnan(MAT_INmean)); % single value

OUT = [mean(cnt_data_ndays)... % mean coverage in n images with spatial binning gap filling
        mean(cnt_data_ndays./npixelsMAR_in)... % mean fractional coverage prior to temporal binning
        cnt_data_ndays_mean...  % mean coverage after temporal binning
        cnt_data_ndays_mean/max(npixelsMAR_in)]; % mean fractional coverage after temporal binning. Need to use MAX number of marine pixels to avoid fractions >1!!!