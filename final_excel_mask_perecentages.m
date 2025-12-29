clc; clear; close all;

%% ================== User Inputs ==================
fat_frac_path   = 'ff_0.nii';                 % fat fraction map (0–1)
short_frac_path = 'Wshort_map_sun4.nii.gz';   % short water fraction map (0–100)

mask_files = {
    'leftQuadMask.nii';
    'rightQuadMask.nii';
    'leftHamsMask.nii';
    'rightHamsMask.nii';
    'VL-R.nii';
    'VL-L.nii';
    'VI-R.nii';
    'VI-L.nii';
    'VM-L.nii';
    'VM-R.nii';
    'RF-R.nii';
    'RF-L.nii';
    'SM-R.nii';
    'SM-L.nii';
    'ST-R.nii';
    'ST-L.nii';
    'BF-R.nii';
    'BF-L.nii';
    'BFS-R.nii';
    'BFS-L.nii'
};

out_excel = 'muscle_quantification.xlsx';

%% ================== Load Maps ==================
fat_frac   = niftiread(fat_frac_path);      % 0–1
short_frac = niftiread(short_frac_path);    % 0–100

muscle_frac = 1 - fat_frac;                 % water fraction
long_frac   = 100 - short_frac;             % long water (%)

%% ================== Prepare Output Table ==================
results = table( ...
    'Size',[numel(mask_files) 8], ...
    'VariableTypes',{'string','double','double','double','double','double','double','double'}, ...
    'VariableNames',{ ...
        'MaskName', ...
        'Fat_pct', ...
        'Water_pct', ...
        'Water_Long_pct', ...
        'Water_Short_pct', ...
        'Final_FA_pct', ...
        'Final_Water_Long_pct', ...
        'Final_Water_Short_pct'} ...
);

%% ================== Loop Over Masks ==================
for i = 1:numel(mask_files)

    maskName = mask_files{i};
    mask     = niftiread(maskName) > 0;

    % ----- Fat / Water -----
    fat_vals    = fat_frac(mask & ~isnan(fat_frac));
    water_vals = muscle_frac(mask & ~isnan(muscle_frac));

    fat_pct   = mean(fat_vals)   * 100;
    water_pct = mean(water_vals) * 100;

    % ----- Short / Long Water -----
    short_vals = short_frac(mask & ~isnan(short_frac));
    long_vals  = long_frac(mask  & ~isnan(long_frac));

    short_pct = mean(short_vals);
    long_pct  = mean(long_vals);

    % ----- Combined Metrics (your exact logic) -----
    final_fa_pct          = fat_pct;
    final_water_long_pct  = (water_pct * long_pct) / 100;
    final_water_short_pct = 100 - final_fa_pct - final_water_long_pct;

    % ----- Store -----
    results.MaskName(i)               = string(erase(maskName,'.nii'));
    results.Fat_pct(i)                = fat_pct;
    results.Water_pct(i)              = water_pct;
    results.Water_Long_pct(i)         = long_pct;
    results.Water_Short_pct(i)        = short_pct;
    results.Final_FA_pct(i)           = final_fa_pct;
    results.Final_Water_Long_pct(i)   = final_water_long_pct;
    results.Final_Water_Short_pct(i)  = final_water_short_pct;
end

%% ================== Write Excel ==================
writetable(results, out_excel);
