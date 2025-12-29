clc; clear; close all;

%% -------------------- User inputs ------------------------
fat_frac_path = 'ff_0.nii';  % short fraction map
%long_mag_path   = 'long_T2star_map.nii.gz';             % long magnitude map (or can be 1-short if fractional)
mask_path       = 'VL_dseg.nii';              % mask of interest

%% -------------------- Load data --------------------------
fat_frac = niftiread(fat_frac_path);  % assumed 0-1
mask       = niftiread(mask_path);

% If you need long fraction explicitly
muscle_frac = 1 - fat_frac;

%% -------------------- Mask the values --------------------
%short_masked = short_frac(mask>0);
%long_masked  = long_frac(mask>0);

fat_masked = fat_frac(mask>0 & ~isnan(fat_frac));
muscle_masked  = muscle_frac(mask>0 & ~isnan(muscle_frac));

%% -------------------- Compute percentages ----------------
fat_pct = mean(fat_masked)*100 ;
muscle_pct  = mean(muscle_masked)*100 ;
fprintf('----Fat/Muscle quantification based on ff map----\n');
fprintf('Within the VL mask:\n');
fprintf('Fat fraction: %.2f %%\n', fat_pct);
fprintf('Muscle fraction:  %.2f %%\n', muscle_pct);



%% -------------------- wf maps ----------------

short_frac_path = 'Wshort_map_mon1.nii.gz';  % short fraction map
%long_mag_path   = 'long_T2star_map.nii.gz';             % long magnitude map (or can be 1-short if fractional)
mask_path2       = 'VL_dseg.nii';              % mask of interest

%% -------------------- Load data --------------------------
short_frac = niftiread(short_frac_path);  % assumed 0-1
mask       = niftiread(mask_path2);

% If you need long fraction explicitly
long_frac = 100 - short_frac;

%% -------------------- Mask the values --------------------
%short_masked = short_frac(mask>0);
%long_masked  = long_frac(mask>0);

short_masked = short_frac(mask>0 & ~isnan(short_frac));
long_masked  = long_frac(mask>0 & ~isnan(long_frac));

%% -------------------- Compute percentages ----------------
short_pct = mean(short_masked) ;
long_pct  = mean(long_masked) ;

fprintf('\n');
fprintf('----Long/Short quantification based on Wshort map----\n');
fprintf('Within the VL mask:\n');
fprintf('Short fraction: %.2f %%\n', short_pct);
fprintf('Long fraction:  %.2f %%\n', long_pct);

%% -------------------- All maps combined ----------------

fprintf('\n');
fprintf('----All maps put together----\n');

fprintf('Fat fraction: %.2f %%\n', fat_pct);

total_Muscle_fraction = (muscle_pct*long_pct)/100;
total_short_fraction = 100-fat_pct-total_Muscle_fraction;
fprintf('Muscle/Long fraction: %.2f %%\n', total_Muscle_fraction);
fprintf('Short fraction: %.2f %%\n', total_short_fraction);