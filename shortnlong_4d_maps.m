clc; clear;

%% ===================== LOAD REFERENCE NIFTI INFO =====================
info = niftiinfo('UTE_mag_orignal_0.nii');

%% ===================== LOAD PARAMETER MAPS =====================
W        = single(niftiread('Water_long_randy_original_mon1.nii.gz'));
F        = single(niftiread('water_short_randy_orignal_mon1.nii.gz'));

% R2* MAPS ARE IN ms^-1
R2water_ms = single(niftiread('R2water_long_map_mon1.nii.gz'));   
R2fat_ms   = single(niftiread('R2water_short_map_mon1.nii.gz'));  

%% ===================== CONVERT R2* TO s^-1 =====================
R2water = R2water_ms * 1000;   % ms^-1 → s^-1
R2fat   = R2fat_ms   * 1000;   % ms^-1 → s^-1

%% ===================== LOAD & CONVERT TEs =====================
t_ms = load('TEs.txt');        % ms
t    = t_ms(:).' / 1000;       % seconds
Nt   = numel(t);

%% ===================== SANITY CHECKS =====================
assert(isequal(size(W), size(F), size(R2water), size(R2fat)), ...
       'All parameter maps must have identical dimensions.');
assert(all(t > 0), 'Echo times must be positive.');

%% ===================== PREALLOCATE =====================
[Nx, Ny, Nz] = size(W);

longComp  = zeros(Nx, Ny, Nz, Nt, 'single');
shortComp = zeros(Nx, Ny, Nz, Nt, 'single');
S         = zeros(Nx, Ny, Nz, Nt, 'single');

%% ===================== SIGNAL GENERATION =====================
for k = 1:Nt
    longComp(:,:,:,k)  = W .* exp(-R2water .* t(k));
    shortComp(:,:,:,k) = F .* exp(-R2fat   .* t(k));
    S(:,:,:,k)         = longComp(:,:,:,k) + shortComp(:,:,:,k);
end

%% ===================== UPDATE NIFTI HEADER =====================
info.ImageSize     = size(S);
info.Datatype      = 'single';
info.BitsPerPixel  = 32;

% Expand pixel dimensions to 4D
if numel(info.PixelDimensions) == 3
    info.PixelDimensions(4) = mean(diff(t));   % TE spacing (seconds)
end

info.raw.dim(1)    = 4;
info.raw.dim(2:5)  = info.ImageSize;

%% ===================== SAVE FILES =====================
niftiwrite(longComp,  'signal_long_4D.nii.gz',  info, 'Compressed', true);
niftiwrite(shortComp, 'signal_short_4D.nii.gz', info, 'Compressed', true);
niftiwrite(S,         'signal_total_4D.nii.gz', info, 'Compressed', true);

disp('4D signals saved with correct units and NIfTI header.');
