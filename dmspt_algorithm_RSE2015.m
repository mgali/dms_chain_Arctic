% Calculate DMSPt concentration using algorithm described in
% Galí et al. 2015 Remote Sensing of Environment 171:171-184
% Accepts vector and N-dimension matrix inputs.
% Consistency in variable size is not tested, Matlab will complain.

% UNITS: chl_in in [mg m-3], zeu_in in [m], mld_in in [m],
% sst_in in [degree celsius], pic_in in ]mol m-3]
% dmspt_out in [micromol m-3]

function [dmspt_out,flags_out,flag_key] = dmspt_algorithm_RSE2015(chl_in,zeu_in,mld_in,sst_in,pic_in)

if nargin == 5
    
    % Define algorithm parameters
    ps = [1.70 1.14 0.44 0.063 -0.0024]; % stratified case params
    pm = [1.74 0.81 0.60]; % mixed case params
    pc = [-1.05 -3.19 -0.78]; % [cocco bloom & chl absent] params
    
    % Restrict chlorophyll range
    chl_in(chl_in < 0.04) = 0.04;
    chl_in(chl_in > 60) = 60;
    
    % Assign sub-algorithms (create algo flags)
    flag_key = {'1_strat' '2_mix' '3_cocc'};
    flags_out = nan(size(chl_in));
    flags_out((zeu_in./mld_in >= 1 | pic_in >= 0.0015) & ~isnan(chl_in) & ~isnan(sst_in)) = 1;
    flags_out(flags_out~=1 & zeu_in./mld_in < 1 & ~isnan(chl_in)) = 2;
    flags_out(pic_in >= 0.0015 & isnan(chl_in)) = 3;
    
    % Log-transform
    chl_in = real(log10(chl_in));
    pic_in = real(log10(pic_in));
    
    % Calculate DMSPt
    dmspt_out = nan(size(chl_in));
    dmspt_out(flags_out == 1) = ps(1) + ps(2)*chl_in(flags_out == 1) + ps(3)*(chl_in(flags_out == 1).^2)...
        + ps(4)*sst_in(flags_out == 1) + ps(5)*(sst_in(flags_out == 1).^2);
    dmspt_out(flags_out == 2) = pm(1) + pm(2)*chl_in(flags_out == 2)...
        + pm(3)*(real(log10(zeu_in(flags_out == 2)./mld_in(flags_out == 2))));
    dmspt_out(flags_out == 3) = pc(1) + pc(2)*pic_in(flags_out == 3) + pc(3)*(pic_in(flags_out == 3).^2);
    dmspt_out = 10.^dmspt_out;
    
else
    error('Function requires 5 arguments, but %i arguments were provided',nargin)
end