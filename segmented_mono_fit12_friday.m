clc; clear; close all;

%% ===================== User Inputs ==========================
mag_path   = 'recon_integrated_norm_complex_v4_mag.nii.gz';  % 4D magnitude
TEs_path   = 'TEs.txt';                                      % echo times (ms)
mask_path  = 'VL_dseg.nii';                             % for ROI noise estimation
out_T2s    = 'T2star_map_ricianLL_filteredTEs_gt3ms_mon1.nii.gz';
out_S0     = 'S0_map_ricianLL_filteredTEs_gt3ms_mon1.nii.gz';
out_R2s    = 'R2star_map_ricianLL_filteredTEs_gt3ms_mon1.nii.gz';
%% ===================== Load Input Data ======================
TEs = double(load(TEs_path));     % ms
TEs = TEs(:);

mag_full = double(niftiread(mag_path));  % (x,y,z,ne)
info     = niftiinfo(mag_path);
info2    = niftiinfo('fat_mag_0.nii');

[nx, ny, nz, ne_full] = size(mag_full);

%% ===================== Filter Echoes > 3 ms ======================
keep = TEs > 3;         % only use echos greater than 3 ms
TEs_filt = TEs(keep);

mag = mag_full(:,:,:,keep);   % filtered magnitude (aligned with TEs_filt)
ne  = numel(TEs_filt);

fprintf('Using %d echoes (TE > 3 ms):\n', ne);
disp(TEs_filt');

%% ================= Load mask for ROI noise estimation ==================
if exist(mask_path, 'file')
    roi.mask = logical(niftiread(mask_path));
else
    roi.mask = true(nx, ny, nz);
end

%% ===================== Estimate sigma using FitImageSigma ==================
imData.images = mag;           % filtered 4D magnitude data
imData.TE = TEs_filt / 1000;   % convert to seconds
imData.FieldStrength = 3;      % Tesla

disp("Estimating sigma using FitImageSigma...")
[maps, sigmaFromRoi, sigmaFromFit] = FittingProcessorRician_3Comp.FitImageSigma(imData, roi);

fprintf('Sigma from ROI: %.4f\n', sigmaFromRoi);
fprintf('Sigma from voxel-wise fitting (corrected): %.4f\n', sigmaFromFit);

sigma = sigmaFromFit;

%% ================= Allocate Output Maps ======================
T2s_map = zeros(nx, ny, nz);
S0_map  = zeros(nx, ny, nz);
R2s_map = zeros(nx, ny, nz);

%% ================= Fitting Bounds ===========================
%lb = [0, 0.5];      % S0 >= 0, T2* >= 0.5 ms
%ub = [3, 200];      % S0 <= 3, T2* <= 200 ms
%ub = [2, 200];      % S0 <= 3, T2* <= 200 ms

%% ================= fmincon Options ==========================
opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ... 
    'Display', 'off', ...
    'MaxFunctionEvaluations', 1000, ...
    'OptimalityTolerance', 1e-8);

 % opts = optimoptions('fmincon', ...
 %                'Algorithm', 'interior-point', ...
 %                'InitBarrierParam', 100000, ...
 %                'ScaleProblem', true, ...
 %                'FiniteDifferenceType', 'central', ...
 %                'Display', 'off');

%% ===================== Parallel Rician LL Fitting ==================
disp("Fitting voxel-wise using filtered TEs (TE > 3 ms)...");
tic

if isempty(gcp('nocreate'))
    parpool;
end

sliceProgress = waitbar(0,'Fitting slices...');

parfor ix = 1:nx
    T2s_slice = zeros(1, ny, nz);
    S0_slice  = zeros(1, ny, nz);

    for iy = 1:ny
        for iz = 1:nz

            S = squeeze(mag(ix,iy,iz,:));  % filtered echoes only

            if all(S==0) || max(S) < 0.001
                continue;
            end

            % Initial guess

            S0_init = max(S);
            lb = [0, 0.5];                % S0 >= 0, T2* >= 0.5 ms
            ub = [3*S0_init, 200];           % S0 <= 3 × observed signal
            %S0_init = 0.2;
            %T2_init = 1;   % ms
            T2_init = 1;   % ms
            p0 = [S0_init, T2_init];

            % Rician log-likelihood
            voxelObj = @(p) -RicianLogLik(S, p(1)*exp(-TEs_filt./p(2)), sigma);

            voxelObjSafe = @(p) max(voxelObj(p), 1e-12);

            try
                p = fmincon( voxelObjSafe, p0, [],[],[],[], lb, ub, [], opts);
                S0_slice(1,iy,iz)  = p(1);
                T2s_slice(1,iy,iz) = p(2);
            catch
            end
        end
    end

    T2s_map(ix,:,:) = T2s_slice;
    S0_map(ix,:,:)  = S0_slice;

    if mod(ix,5)==0
        waitbar(ix/nx, sliceProgress, sprintf('Fitting slice %d / %d', ix, nx));
    end
end

close(sliceProgress)
toc

R2s_map = zeros(size(T2s_map));
valid = T2s_map > 0;
R2s_map(valid) = 1 ./ T2s_map(valid);   % ms^-1
%% ===================== Save NIfTI Output =====================
info2.Datatype = 'double';
niftiwrite(T2s_map, out_T2s, info2);
niftiwrite(S0_map,  out_S0, info2);
niftiwrite(R2s_map, out_R2s, info2);

disp("DONE! Saved filtered-TE Rician LL T2* and S0 maps.");

%% ===================== Rician Log-Likelihood Function =====================
function [loglik, logliks] = RicianLogLik(measurements, predictions, sigma)
    measurements = max(measurements, eps);
    sigmaSquared = sigma.^2;
    sumsqsc = (measurements.^2 + predictions.^2) ./ (2*sigmaSquared);
    scp = measurements.*predictions ./ sigmaSquared;
    lb0 = logbesseli0(scp);
    logliks = log(measurements) - log(sigmaSquared) - sumsqsc + lb0;
    loglik = sum(logliks);
end

function lb0 = logbesseli0(x)
    lb0 = zeros(size(x));
    exact = x < 700;
    approx = ~exact;
    lb0(exact) = log(besseli(0,x(exact)));
    lb0(approx) = x(approx) - 0.5*log(2*pi*x(approx));
end
