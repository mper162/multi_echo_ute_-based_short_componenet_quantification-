clc; clear; close all;

%% ================== User Inputs ==================
fat_frac_path    = 'ff_0.nii';
short_frac_path  = 'Wshort_map_mon1.nii.gz';
orig_ute_path    = 'UTE_mag_orignal_0.nii';

t2_short_path = 'T2water_short_map_mon1.nii.gz';
t2_long_path  = 'T2water_long_map_mon1.nii.gz';

mask_files = {
    'leftQuadMask.nii';
    'rightQuadMask.nii';
    'leftHamsMask.nii';
    'rightHamsMask.nii';
    'VL-R.nii'; 'VL-L.nii';
    'VI-R.nii'; 'VI-L.nii';
    'VM-L.nii'; 'VM-R.nii';
    'RF-R.nii'; 'RF-L.nii';
    'SM-R.nii'; 'SM-L.nii';
    'ST-R.nii'; 'ST-L.nii';
    'BF-R.nii'; 'BF-L.nii';
    'BFS-R.nii'; 'BFS-L.nii'
};

out_excel = 'muscle_quantification_30.xlsx';

%% ================== Load Maps (FORCE DOUBLE) ==================
fat_frac     = double(niftiread(fat_frac_path));
short_frac   = double(niftiread(short_frac_path));
t2_short_map = double(niftiread(t2_short_path));
t2_long_map  = double(niftiread(t2_long_path));

muscle_frac = 1 - fat_frac;
long_frac   = 100 - short_frac;

%% ================== Load Original Image Header ==================
ute_info = niftiinfo(orig_ute_path);
voxel_volume_mL = prod(ute_info.PixelDimensions(1:3)) * 1e-3;

%% ================== Output Table ==================
results = table( ...
    'Size',[numel(mask_files) 11], ...
    'VariableTypes',{ ...
        'string','double','double','double','double', ...
        'double','double','double','double','double','double'}, ...
    'VariableNames',{ ...
        'MaskName','Volume_mL','Fat_pct','Water_pct', ...
        'Water_Long_pct','Water_Short_pct', ...
        'Final_FA_pct','Final_Water_Long_pct','Final_Water_Short_pct', ...
        'Water_Long_T2_ms','Water_Short_T2_ms'} ...
);

%% ================== Loop Over Masks ==================
for i = 1:numel(mask_files)

    maskName = mask_files{i};
    mask = niftiread(maskName) > 0;

    % ---------- Volume ----------
    n_voxels = nnz(mask);
    volume_mL = n_voxels * voxel_volume_mL;

    % ---------- Fat / Water ----------
    fat_vals   = fat_frac(mask);
    water_vals = muscle_frac(mask);

    fat_vals   = fat_vals(isfinite(fat_vals) & fat_vals >= 0 & fat_vals <= 1);
    water_vals = water_vals(isfinite(water_vals) & water_vals >= 0 & water_vals <= 1);

    fat_pct   = mean(fat_vals) * 100;
    water_pct = mean(water_vals) * 100;

    % ---------- Short / Long Fractions ----------
    short_vals = short_frac(mask);
    long_vals  = long_frac(mask);

    short_vals = short_vals(isfinite(short_vals) & short_vals > 0 & short_vals < 100);
    long_vals  = long_vals(isfinite(long_vals)  & long_vals > 0 & long_vals < 100);

    short_pct = mean(short_vals);
    long_pct  = mean(long_vals);

    % ---------- T2* VALUES (CRITICAL FIX) ----------
    t2s = t2_short_map(mask);
    t2l = t2_long_map(mask);

    % Remove Inf, NaN, zero, and 65535
    valid_short = isfinite(t2s) & t2s > 0 & t2s < 100 & t2s ~= 65535;
    valid_long  = isfinite(t2l) & t2l > 0 & t2l < 200 & t2l ~= 65535;

    mean_t2_short = mean(t2s(valid_short));
    mean_t2_long  = mean(t2l(valid_long));

    % ---------- Combined Metrics ----------
    final_fa_pct          = fat_pct;
    final_water_long_pct  = (water_pct * long_pct) / 100;
    final_water_short_pct = 100 - final_fa_pct - final_water_long_pct;

    % ---------- Store ----------
    results.MaskName(i)              = string(erase(maskName,'.nii'));
    results.Volume_mL(i)             = volume_mL;
    results.Fat_pct(i)               = fat_pct;
    results.Water_pct(i)             = water_pct;
    results.Water_Long_pct(i)        = long_pct;
    results.Water_Short_pct(i)       = short_pct;
    results.Final_FA_pct(i)          = final_fa_pct;
    results.Final_Water_Long_pct(i)  = final_water_long_pct;
    results.Final_Water_Short_pct(i) = final_water_short_pct;
    results.Water_Long_T2_ms(i)      = mean_t2_long;
    results.Water_Short_T2_ms(i)     = mean_t2_short;

end

%% ================== Write Excel ==================
writetable(results, out_excel);
