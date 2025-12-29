%% full_recon_subtract_fat_normalized_complex_v3_autoMask.m
% Complex UTE fat subtraction using:
%  - 3D fat magnitude + 3D fat phase
%  - Dixon multi-peak fat model
%  - B0 map
%  - Automatic fat mask via thresholding
%  - Complex LS scaling per echo
%
% OUTPUT:
%   recon_integrated_norm_complex_v3_* (magnitude + phase)

clc; clear; close all;

fprintf('\n===== UTE + FAT COMPLEX SUBTRACTION (AUTO FAT MASK) =====\n');

%% -------------------- User inputs ------------------------
ute_mag_path      = 'UTE_mag_orignal_0.nii';      % 4D UTE magnitude
ute_phase_path    = 'UTE_phase_original_0.nii';  % 4D UTE phase
fat_mag_3d_path   = 'fat_mag_0.nii';              % 3D fat magnitude
%fat_phase_3d_path = 'fat_phase_0.nii';            % 3D fat phase
b0_path           = 'b0_0.nii';                    % 3D B0 map (Hz)
r2star_path       = 'r2s_0.nii';                % optional R2* map (s^-1)
TEs_path          = 'TEs.txt';                     % echo times (ms or s)

TE_ref_index      = 1;                             % reference echo
clipPct           = 99.9;                          % fat magnitude clipping
out_prefix        = 'recon_integrated_norm_complex_v4';

% Auto fat mask parameters
fat_thr_pct = 85;      % percentile for fat-dominant voxels
ute_thr_pct = 20;      % background suppression

tesla = 3.0;

%% -------------------- Load inputs ------------------------
fprintf('Loading data...\n');

ute_info = niftiinfo(ute_mag_path);
ute_mag4 = double(niftiread(ute_mag_path));
ute_ph4  = double(niftiread(ute_phase_path));

fat_mag3 = double(niftiread(fat_mag_3d_path));
%fat_ph3  = double(niftiread(fat_phase_3d_path));
b0map    = double(niftiread(b0_path));

TEs = load(TEs_path);
if max(TEs) > 1
    TEs = TEs / 1000;  % ms → s
end

[xdim, ydim, zdim, ne] = size(ute_mag4);
assert(ne == numel(TEs), 'Echo count mismatch');

%% -------------------- Phase normalization ------------------------
ute_ph4 = normalize_phase(ute_ph4);
%fat_ph3 = normalize_phase(fat_ph3);


%% -------------------- Build complex UTE ------------------------
signal_total = ute_mag4 .* exp(1i * ute_ph4);
fat_ph3 = angle(signal_total(:,:,:,TE_ref_index));
clear ute_mag4 ute_ph4

%% -------------------- Build fat prior at reference TE ------------------------
fat_complex_ref = fat_mag3 .* exp(1i * fat_ph3);

% Clip extreme fat magnitudes
fclip = prctile(abs(fat_complex_ref(:)), clipPct);
fat_prior_mag = min(abs(fat_complex_ref), fclip);

%% -------------------- Dixon multi-peak fat evolution ------------------------
fprintf('Building Dixon-consistent fat signal...\n');

dTEs = TEs - TEs(TE_ref_index);
fat_peaks = DixonMultiPeakFat(dTEs, tesla);   % [10 x NE]

Fprior_4d = complex(zeros(xdim, ydim, zdim, ne));

for e = 1:ne
    fat_sum = sum(fat_peaks(:,e));
    b0_rot  = exp(1i * 2*pi * b0map * dTEs(e));

    Fprior_4d(:,:,:,e) = ...
        fat_prior_mag .* ...
        exp(1i * fat_ph3) .* ...
        fat_sum .* ...
        b0_rot;
end

%% -------------------- Automatic fat mask ------------------------
fprintf('Estimating fat-dominant voxels...\n');

fat_mag_ref = fat_prior_mag;
ute_mag_ref = abs(signal_total(:,:,:,TE_ref_index));

fat_thr = prctile(fat_mag_ref(:), fat_thr_pct);
ute_thr = prctile(ute_mag_ref(:), ute_thr_pct);

fat_mask = fat_mag_ref > fat_thr & ute_mag_ref > ute_thr;

fprintf('Auto fat mask voxels: %d (%.2f%%)\n', ...
    nnz(fat_mask), 100 * nnz(fat_mask) / numel(fat_mask));

%% -------------------- Complex LS fat subtraction ------------------------
fprintf('Performing complex LS subtraction...\n');

recon_complex_final = complex(zeros(xdim, ydim, zdim, ne));
fat_scaled_final   = complex(zeros(xdim, ydim, zdim, ne));

for e = 1:ne
    S_ute = signal_total(:,:,:,e);
    S_fat = Fprior_4d(:,:,:,e);

    A = S_fat(fat_mask);
    B = S_ute(fat_mask);

    scale = (A(:)' * B(:)) / (A(:)' * A(:) + eps);

    S_fat_scaled = scale * S_fat;

    recon_complex_final(:,:,:,e) = S_ute - S_fat_scaled;
    fat_scaled_final(:,:,:,e)   = S_fat_scaled;
end

%% -------------------- Save outputs ------------------------
fprintf('Saving NIfTI outputs...\n');

info4d = ute_info;
info4d.Datatype = 'single';

niftiwrite(single(abs(recon_complex_final)), ...
    [out_prefix '_mag.nii'], info4d, 'Compressed', true);

niftiwrite(single(angle(recon_complex_final)), ...
    [out_prefix '_phase.nii'], info4d, 'Compressed', true);

niftiwrite(single(abs(fat_scaled_final)), [out_prefix '_fatmag.nii'], info4d, 'Compressed', true);
niftiwrite(single(angle(fat_scaled_final)), [out_prefix '_fatphase.nii'], info4d, 'Compressed', true);

fprintf('Done. Output prefix: %s\n', out_prefix);

%% ==================== FUNCTIONS ====================

function ph = normalize_phase(ph_raw)
    p99 = prctile(ph_raw(:),99);
    p1  = prctile(ph_raw(:),1);
    rngv = p99 - p1;

    if all(mod(ph_raw(:),1)==0) && max(ph_raw(:)) >= 1024
        ph = (ph_raw / double(max(ph_raw(:)))) * 2*pi - pi;
    elseif rngv > 2*pi && rngv < 360
        ph = ph_raw * (pi/180);
    else
        ph = ph_raw;
    end
end

function fat_rel = DixonMultiPeakFat(TEs, tesla)
    gyro = 42.58e6;
    larmor = tesla * gyro;

    FatPPM = [-3.81 -3.4 -3.12 -2.67 -2.45 -1.94 -0.63 -0.4 0.52 0.62];
    FatAmp = [0.089 0.577 0.059 0.093 0.059 0.013 0.02 0.02 0.01 0.059];
    FatR2  = [130 135 85 95 100 90 100 100 85 85];

    FatFreq = FatPPM * larmor * 1e-6;
    FatW = 2*pi*FatFreq;

    TEs = TEs(:).';
    nE = numel(TEs);
    nP = numel(FatAmp);

    fat_rel = complex(zeros(nP, nE));

    for k = 1:nP
        fat_rel(k,:) = ...
            FatAmp(k) .* ...
            exp(1i * FatW(k) * TEs) .* ...
            exp(-FatR2(k) * abs(TEs));
    end
end
