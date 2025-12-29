classdef FittingProcessorRician_3Comp
    % FittingProcessorRician_3Comp
    %   This class implements multi‐step image fitting routines for water–fat
    %   separation using a three‐component model. In this model, the signal is
    %   assumed to be a sum of fat, water_long and water_short components – each with
    %   its own R2 decay (R2fat, R2water_long, and R2water_short). The initialization
    %   for the two‐point search forces an equal split between water_long and water_short,
    %   but once the fitting is done the fitted water components are used as obtained.
    %   For sigma estimation the usual two‐component model is used.
    
    methods (Static)
        function [maps, sigmaFromRoi, sigmaFromFit] = FitImageSigma(imData, roi)
            % unchanged from your original code
            TE = 1000 * imData.TE;
            TE = reshape(TE, [], 1);
            FieldStrength = imData.FieldStrength;
            
            mask = roi.mask;
            data = imData.images;
            [nx, ny, nz, nTE] = size(data);
            
            % Preallocate maps as column vectors
            FFrician = zeros(nx*ny*nz, 1);
            R2water_rician = zeros(nx*ny*nz, 1);
            R2fat_rician   = zeros(nx*ny*nz, 1);
            sigmaEstimates = zeros(nx*ny*nz, 1);
            s0Estimates    = zeros(nx*ny*nz, 1);
            
            mask_r = mask(:);
            data_r = reshape(data, nx*ny*nz, nTE);
            RoiIdxs = mask_r > 0 & all(data_r, 2);
            data_roi = data_r(RoiIdxs, :);
            nb_voxels = size(data_roi, 1);
            
            % Estimate sigma from ROI using standard deviation across echoes.
            stdMaskedIm = std(data_roi, [], 1);
            [~, TEindex] = min(abs(imData.TE - 0.0023));
            sigmaFromRoi = stdMaskedIm(TEindex);
            
            % For sigma estimation we use option 2 (2-component model)
            algoparamsArray = FittingProcessorRician_3Comp.setAlgoparamsArray(data_roi, TE, sigmaFromRoi, 2);
            
            FFrician_roi = zeros(nb_voxels, 1);
            R2water_roi   = zeros(nb_voxels, 1);
            R2fat_roi     = zeros(nb_voxels, 1);
            sigmaEstimates_roi = zeros(nb_voxels, 1);
            s0Estimates_roi    = zeros(nb_voxels, 1);
            
            disp('Estimating sigma from ROI (2-component model)...');
            tstart = tic;
            parfor i = 1:nb_voxels
                ydata = data_roi(i,:);
                algoparams = algoparamsArray(i);
                outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting_WithSigma(TE, FieldStrength, ydata, algoparams.sigEst, algoparams);
                % The 2-component model returns parameters: [F; W; R2; sigma]
                FFrician_roi(i) = outparams.F / (outparams.F + outparams.W);
                R2water_roi(i) = outparams.R2water;
                R2fat_roi(i)   = outparams.R2fat;
                sigmaEstimates_roi(i) = outparams.sig;
                s0Estimates_roi(i) = outparams.F + outparams.W;
            end
            tstop = toc(tstart);
            disp("Sigma estimated in: " + tstop + "s");
            
            FFrician(RoiIdxs,:) = FFrician_roi;
            R2water_rician(RoiIdxs,:) = R2water_roi;
            R2fat_rician(RoiIdxs,:) = R2fat_roi;
            sigmaEstimates(RoiIdxs,:) = sigmaEstimates_roi;
            s0Estimates(RoiIdxs,:) = s0Estimates_roi;
            
            maps.FFrician = reshape(FFrician, nx, ny, nz);
            maps.R2water   = reshape(R2water_rician, nx, ny, nz);
            maps.R2fat     = reshape(R2fat_rician, nx, ny, nz);
            maps.sigma     = reshape(sigmaEstimates, nx, ny, nz);
            maps.S0        = reshape(s0Estimates, nx, ny, nz);
            maps.SNR       = s0Estimates ./ sigmaEstimates;
            
            validSigma = maps.sigma(maps.sigma > 0);
            sigmaFromFit = mean(validSigma);
            correctionData = 1.163;
            sigmaFromFit = sigmaFromFit * correctionData;
        end


%%%%%%%%Randy first part%%%%%%%%%%%%%%%%%%%%

        function maps = FitImageSigma_randy(imData, sigma)
            % unchanged from your original code
            TE = 1000 * imData.TE;
            TE = reshape(TE, [], 1);
            FieldStrength = imData.FieldStrength;
            data = imData.images;
            [nx, ny, nz, nTE] = size(data);

            %data_long = imData.Water_long;
            %data_r2_long  = imData.R2water_long;

           % [nx_long, ny_long, nz_long] = size(data_long);
           % [nx_r2_long, ny_r2_long, nz_r2_long] = size(data_r2_long);

            if isfield(imData, 'mask')
                mask = imData.mask;
            else
                Img2norm = data(:,:,:,1);
                threshold = prctile(Img2norm(:), 85);
                mask = Img2norm > threshold;
            end

            tmp = squeeze(sum(sum(mask,1),2));
            slice_idxs = find(tmp > 0);
            
            % Preallocate maps as column vectors
            FFrician = zeros(nx*ny*nz, 1);
            R2water_rician = zeros(nx*ny*nz, 1);
            R2fat_rician   = zeros(nx*ny*nz, 1);
            sigmaEstimates = zeros(nx*ny*nz, 1);
            s0Estimates    = zeros(nx*ny*nz, 1);
            Wl_est = zeros(nx*ny*nz,1);
            Ws_est = zeros(nx*ny*nz,1);

            mask_r = mask(:);
            data_r = reshape(data, nx*ny*nz, nTE);
            RoiIdxs = mask_r > 0 & all(data_r, 2);
            data_roi = data_r(RoiIdxs, :);
            nb_voxels = size(data_roi, 1);

            %%%%% water_Long fixed
            % if isfield(imData, 'mask')
            %     mask_long = imData.mask;
            % else
            %     Img2norm_long =  data_long;
            %     threshold_long = prctile(Img2norm_long(:), 85);
            %     mask_long = Img2norm_long > threshold_long;
            % end
            % 
            % tmp_long = squeeze(sum(sum(mask_long,1),2));  % sum over x and y
            % slice_idxs_long = find(tmp_long > 0);          % indices of slices with signal
            % 
            % map_long_r = data_long(:);  % flatten 3D map
            % mask_long_r = mask_long(:);               % flatten mask
            % 
            % %data_long_r = reshape(data_long, nx_long*ny_long*nz_long, nTE);
            % % RoiIdxs_long = mask_long_r > 0 & all(data_long_r, 2);
            % % data_long_roi = data_long_r(RoiIdxs_long, :);
            % % nb_voxels_long = size(data_long_roi, 1);
            % 
            % RoiIdxs_long = mask_long_r > 0 & ~isnan(map_long_r) & map_long_r ~= 0;  % only positive & non-NaN voxels
            % map_roi_long = map_long_r(RoiIdxs_long);                            % extract ROI voxels
            % nb_voxels_long = length(map_roi_long);   
            % 
            %%% end of long fixed

            %%%%% water_R2 fixed************************
            % if isfield(imData, 'mask')
            %     mask_r2_long = imData.mask;
            % else
            %     Img2norm_r2_long =  data_r2_long;
            %     threshold_r2_long = prctile(Img2norm_r2_long(:), 85);
            %     mask_r2_long = Img2norm_r2_long > threshold_r2_long;
            % end
            % 
            % tmp_r2_long = squeeze(sum(sum(mask_r2_long,1),2));  % sum over x and y
            % slice_idxs_r2_long = find(tmp_r2_long > 0);          % indices of slices with signal
            % 
            % map_r2_long_r = data_r2_long(:);  % flatten 3D map
            % mask_r2_long_r = mask_r2_long(:);               % flatten mask
            % 
            % %data_long_r = reshape(data_long, nx_long*ny_long*nz_long, nTE);
            % % RoiIdxs_long = mask_long_r > 0 & all(data_long_r, 2);
            % % data_long_roi = data_long_r(RoiIdxs_long, :);
            % % nb_voxels_long = size(data_long_roi, 1);
            % 
            % RoiIdxs_r2_long = mask_r2_long_r > 0 & ~isnan(map_r2_long_r) & map_r2_long_r ~= 0;  % only positive & non-NaN voxels
            % map_roi_r2_long = map_r2_long_r(RoiIdxs_r2_long);                            % extract ROI voxels
            % nb_voxels_r2_long = length(map_roi_r2_long);   
            
            %%% end of R2 fixed***************************


            disp('Setting fit parameters for 3-component model...');
            tstart = tic;
            algoparamsArray = FittingProcessorRician_3Comp.setAlgoparamsArray(data_roi, TE, sigma,4);
            tstop = toc(tstart);
            disp("Parameters set in " + tstop + " seconds");
            
            FFrician_roi = zeros(nb_voxels, 1);
            R2water_roi   = zeros(nb_voxels, 1);
            R2fat_roi     = zeros(nb_voxels, 1);
            sigmaEstimates_roi = zeros(nb_voxels, 1);
            s0Estimates_roi    = zeros(nb_voxels, 1);
            Wl_roi = zeros(nb_voxels,1);
            Ws_roi = zeros(nb_voxels,1);

            D = parallel.pool.DataQueue;
            progressCount = 0;
            hWaitbar = waitbar(0, 'Performing 2-component mapping...');
            afterEach(D, @updateWaitbar);
            function updateWaitbar(~)
                progressCount = progressCount + 1;
                waitbar(progressCount/nb_voxels, hWaitbar, ...
                    sprintf('Progress: %d%%', round(100*progressCount/nb_voxels)));
                drawnow;
            end
            

            tstart = tic;
            parfor i = 1:nb_voxels
                ydata = data_roi(i,:);
                algoparams = algoparamsArray(i);
                outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting_WithSigma_randy(TE, FieldStrength, ydata, algoparams.sigEst, algoparams);
                %outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting_FixLong(TE, FieldStrength, ydata, algoparams.sigEst, algoparams);
                
                % The 2-component model returns parameters: [F; W; R2; sigma]
                FFrician_roi(i) = outparams.F / (outparams.F + outparams.W);
                R2water_roi(i) = outparams.R2water;
                R2fat_roi(i)   = outparams.R2fat;
                %sigmaEstimates_roi(i) = outparams.sig;
                s0Estimates_roi(i) = outparams.F + outparams.W;
                Wl_roi(i) = outparams.W;
                Ws_roi(i) = outparams.F;
                send(D, i);
            end
            tstop = toc(tstart);
            delete(hWaitbar);
            disp("2-component data fitted in " + tstop + " seconds");
            
            FFrician(RoiIdxs,:) = FFrician_roi;
            R2water_rician(RoiIdxs,:) = R2water_roi;
            R2fat_rician(RoiIdxs,:) = R2fat_roi;
            %sigmaEstimates(RoiIdxs,:) = sigmaEstimates_roi;
            s0Estimates(RoiIdxs,:) = s0Estimates_roi;
            Wl_est(RoiIdxs,:) = Wl_roi;
            Ws_est(RoiIdxs,:) = Ws_roi;

            maps.FFrician = reshape(FFrician, nx, ny, nz);
            maps.R2water   = reshape(R2water_rician, nx, ny, nz);
            maps.R2fat     = reshape(R2fat_rician, nx, ny, nz);
            %maps.sigma     = reshape(sigmaEstimates, nx, ny, nz);
            maps.S0        = reshape(s0Estimates, nx, ny, nz);
            %maps.SNR       = s0Estimates ./ sigmaEstimates;
            maps.Water_long = reshape(Wl_est, nx, ny, nz);
            maps.Water_short = reshape(Ws_est, nx, ny, nz);
        end

    %%%%%%% end of Randy First Part%%%%%%%%%%%

    
        function maps = FitImage(imData, sigma)
            % FitImage performs voxel-wise fitting on a given slice using a fixed sigma.
            % This routine estimates the three-component model:
            %   [F; Wl; Ws; R2fat; R2water_long; R2water_short]
            % or a reduced model when long-water values are provided:
            %   free parameters [F; Ws; R2fat; R2water_short] with Wl and R2water_long fixed.
            
            %TE = 1000 * imData.TE;
            if max(imData.TE) < 0.1
                TE = 1000 * imData.TE;  % s → ms
            else
                TE = imData.TE;         % already ms
            end

            TE = reshape(TE, [], 1);
            FieldStrength = imData.FieldStrength;
            data = imData.images;
            [nx, ny, nz, nTE] = size(data);
            
            if isfield(imData, 'mask')
                mask = imData.mask;
            else
                Img2norm = data(:,:,:,1);
                threshold = prctile(Img2norm(:), 85);
                mask = Img2norm > threshold;
            end
            
            tmp = squeeze(sum(sum(mask,1),2));
            slice_idxs = find(tmp > 0);
            
            % Preallocate output vectors.
            Wsrician = zeros(nx*ny*nz,1);
            R2fat_rician = zeros(nx*ny*nz,1);
            R2water_long_rician = zeros(nx*ny*nz,1);
            R2water_short_rician = zeros(nx*ny*nz,1);
            s0rician = zeros(nx*ny*nz,1);
            Wl_est = zeros(nx*ny*nz,1);
            Ws_est = zeros(nx*ny*nz,1);
            F_est = zeros(nx*ny*nz,1);
            FFrician = zeros(nx*ny*nz,1);

            mask_r = mask(:);
            data_r = reshape(data, nx*ny*nz, nTE);
            RoiIdxs = mask_r > 0 & all(data_r > 0, 2);
            data_roi = data_r(RoiIdxs, :);
            nb_voxels = size(data_roi, 1);
            
            disp('Setting fit parameters for 2 or 3-component model...');
            tstart = tic;
            
            % --- Detect whether user provided Fixed long-water maps ---
            hasFixedLong = isfield(imData, 'Water_long') && isfield(imData, 'R2water_long');
            if hasFixedLong
                % maps must be same 3D size as image volumes
                Wl_map = imData.Water_long;
                R2wl_map = imData.R2water_long;
                Wl_map_r = Wl_map(:);
                R2wl_map_r = R2wl_map(:);
                % keep only values for RoiIdxs
                fixedWl_roi = Wl_map_r(RoiIdxs);
                fixedR2wl_roi = R2wl_map_r(RoiIdxs);
                % Build algoparamsArray for reduced-fit (opt==3)
                algoparamsArray = FittingProcessorRician_3Comp.setAlgoparamsArray(data_roi, TE, sigma, 3, fixedWl_roi, fixedR2wl_roi);
            else
                algoparamsArray = FittingProcessorRician_3Comp.setAlgoparamsArray(data_roi, TE, sigma, 1);
            end
            
            tstop = toc(tstart);
            disp("Parameters set in " + tstop + " seconds");
            
            % Preallocate results for each voxel.
            Wsrician_roi = zeros(nb_voxels,1);
            R2fat_roi = zeros(nb_voxels,1);
            R2water_long_roi = zeros(nb_voxels,1);
            R2water_short_roi = zeros(nb_voxels,1);
            s0rician_roi = zeros(nb_voxels,1);
            Wl_roi = zeros(nb_voxels,1);
            Ws_roi = zeros(nb_voxels,1);
            F_roi = zeros(nb_voxels,1);
            FFrician_roi = zeros(nb_voxels,1);
            
            disp('Fitting 3-component data...');
            % Set up a progress bar via DataQueue and waitbar.
            D = parallel.pool.DataQueue;
            progressCount = 0;
            hWaitbar = waitbar(0, 'Performing 3-component mapping...');
            afterEach(D, @updateWaitbar);
            function updateWaitbar(~)
                progressCount = progressCount + 1;
                waitbar(progressCount/nb_voxels, hWaitbar, ...
                    sprintf('Progress: %d%%', round(100*progressCount/nb_voxels)));
                drawnow;
            end
            
            tstart = tic;
            if hasFixedLong
                % Reduced fit: fix Wl and R2water_long per-voxel, fit [F; Ws; R2fat; R2water_short]
                parfor i = 1:nb_voxels
                    ydata = data_roi(i,:);
                    algoparams = algoparamsArray(i);
                    %outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting_FixMuscle(TE, FieldStrength, ydata, algoparams.sigEst, algoparams);
                    outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting_FixLong(TE, FieldStrength, ydata, algoparams.sigEst, algoparams); 
                   % den = outparams.F + outparams.Wl + outparams.Ws + eps;
                    den = outparams.Wl + outparams.Ws + eps;
                    %FFrician_roi(i) = outparams.F / den;
                    %FFrician_roi(i) = outparams.F / (outparams.F + outparams.Wl + outparams.Ws);
                    Wsrician_roi(i) = outparams.Ws / den;
                    %R2fat_roi(i) = outparams.R2fat;
                    R2water_long_roi(i) = outparams.R2water_long;
                    R2water_short_roi(i) = outparams.R2water_short;
                    %s0rician_roi(i) =  outparams.Wl + outparams.Ws + outparams.F;
                    s0rician_roi(i) =  outparams.Wl + outparams.Ws;
                    Wl_roi(i) = outparams.Wl;
                    Ws_roi(i) = outparams.Ws;
                    %F_roi(i) = outparams.F;
                    send(D, i);
                end

            else
                % Full fit (unchanged)
                parfor i = 1:nb_voxels
                    ydata = data_roi(i,:);
                    algoparams = algoparamsArray(i);
                    outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting(TE, FieldStrength, ydata, algoparams.sigEst, algoparams);
                    FFrician_roi(i) = outparams.F / (outparams.F + outparams.Wl + outparams.Ws);
                    R2fat_roi(i) = outparams.R2fat;
                    R2water_long_roi(i) = outparams.R2water_long;
                    R2water_short_roi(i) = outparams.R2water_short;
                    s0rician_roi(i) =  outparams.Wl + outparams.Ws + outparams.F;
                    Wl_roi(i) = outparams.Wl;
                    Ws_roi(i) = outparams.Ws;
                    send(D, i);
                end
            end
            tstop = toc(tstart);
            delete(hWaitbar);
            disp("3-component data fitted in " + tstop + " seconds");
            
            % (No post-fit forcing is done; the water split remains as fitted.)
            
            Wsrician(RoiIdxs,:) = Wsrician_roi;
            %R2fat_rician(RoiIdxs,:) = R2fat_roi;
            R2water_long_rician(RoiIdxs,:) = R2water_long_roi;
            R2water_short_rician(RoiIdxs,:) = R2water_short_roi;
            s0rician(RoiIdxs,:) = s0rician_roi;
            Wl_est(RoiIdxs,:) = Wl_roi;
            Ws_est(RoiIdxs,:) = Ws_roi;
            %F_est(RoiIdxs,:) = F_roi;
            %FFrician(RoiIdxs,:) =  FFrician_roi;

            % Assemble output maps.
            maps.Ws = reshape(100.*Wsrician, nx, ny, nz);
            %maps.R2fat = reshape(R2fat_rician, nx, ny, nz);
            maps.R2water_long = reshape(R2water_long_rician, nx, ny, nz);
            maps.R2water_short = reshape(R2water_short_rician, nx, ny, nz);
            %maps.T2str_fat = 1./reshape(R2fat_rician, nx, ny, nz);
            maps.T2str_water_long = 1./reshape(R2water_long_rician, nx, ny, nz);
            maps.T2str_water_short = 1./reshape(R2water_short_rician, nx, ny, nz);
            maps.S0 = reshape(s0rician, nx, ny, nz);
            %maps.Fat = reshape(s0rician .* FFrician, nx, ny, nz);
            %maps.Water = reshape(s0rician - s0rician .* FFrician, nx, ny, nz);
            maps.Water_long = reshape(Wl_est, nx, ny, nz);
            maps.Water_short = reshape(Ws_est, nx, ny, nz);
            %maps.F = reshape(F_est, nx, ny, nz);
            %maps.Wsfraction = 100 .* maps.Water_short ./( maps.Water_short + maps.Water_long + maps.F);
            maps.Wsfraction = 100 .* maps.Water_short ./( maps.Water_short + maps.Water_long);
            %maps.FFrician = 100 .* maps.F ./( maps.Water_short + maps.Water_long + maps.F);
            
            % Display a representative slice.
            % if numel(slice_idxs) > 1
            %     slice = slice_idxs(round(numel(slice_idxs)/2));
            % else
            %     slice = slice_idxs;
            % end
            
            % figure;
            % subplot(2,5,1); imagesc(maps.FFrician(:,:,slice)); colorbar; title('FF (%)'); clim([0 100]);
            % subplot(2,5,2); imagesc(maps.R2fat(:,:,slice)); colorbar; title('R2_{fat} (s^{-1})'); clim([0 250]);
            % subplot(2,5,3); imagesc(maps.R2water_long(:,:,slice)); colorbar; title('R2_{water\_long} (s^{-1})'); clim([0 200]);
            % subplot(2,5,4); imagesc(maps.R2water_short(:,:,slice)); colorbar; title('R2_{water\_short} (s^{-1})'); clim([0 200]);
            % subplot(2,5,5); imagesc(maps.S0(:,:,slice)); colorbar; title('S0');
            % subplot(2,5,6); imagesc(maps.T2str_fat(:,:,slice)); colorbar; title('T2*_{fat} (ms)'); clim([0 80]);
            % subplot(2,5,7); imagesc(maps.T2str_water_long(:,:,slice)); colorbar; title('T2*_{water\_long} (ms)'); clim([0 80]);
            % subplot(2,5,8); imagesc(maps.T2str_water_short(:,:,slice)); colorbar; title('T2*_{water\_short} (ms)'); clim([0 10]);
            % subplot(2,5,9); imagesc(maps.Wsfraction(:,:,slice)); colorbar; title('Water\_short Fraction'); clim([0 100]);
        end
        
        function maps = MultistepFitImage(imData, ~)
            % MultistepFitImage combines sigma estimation and fixed-sigma fitting.
            % It first estimates sigma using the 2-component model and then performs the
            % 3-component fit (or reduced fit if fixed long-water maps provided).
            
              % [~, ~, sigmaFromFit] = FittingProcessorRician_3Comp.FitImageSigma(imData, roi);
              % disp("Estimated sigma is: " + sigmaFromFit);
             % % Use randy variant for final fit (keeps previous behaviour)
              %mapsFittedSigma = FittingProcessorRician_3Comp.FitImageSigma_randy(imData, sigmaFromFit);
             sigmaFromFit =   68.7407;
             % If the user supplied fixed long-water maps they will be honoured in FitImage
             % (we keep mapsFittedSigma as fallback but call FitImage to produce full maps).
             %try
                 maps = FittingProcessorRician_3Comp.FitImage(imData, sigmaFromFit);
             %catch
                 % If anything fails, return the earlier randy result
                 %maps = mapsFittedSigma;
             %end
        end
        
        function algoparamsArray = setAlgoparamsArray(S, echotimes, sigmaEstimateFromRoi, opt, fixedWl_roi, fixedR2wl_roi)
            % setAlgoparamsArray creates an array of fitting parameter structures for each voxel.
            %
            % For the 3-component model (opt == 1), the parameter vector is:
            %   [F; Wl; Ws; R2fat; R2water_long; R2water_short]
            %
            % For sigma estimation (opt == 2), the standard 2-component model is used:
            %   [F; W; R2water; R2fat; sigma]
            %
            % For reduced-fit when long-water is fixed (opt == 3), the free parameter
            % vector is:
            %   [F; Ws; R2fat; R2water_short]
            %
            % Inputs:
            %   fixedWl_roi, fixedR2wl_roi - required for opt==3 (vectors length = nvoxels)
            
            if nargin < 5
                fixedWl_roi = [];
            end
            if nargin < 6
                fixedR2wl_roi = [];
            end
            
            nvoxels = size(S, 1);
            solver = 'fmincon';
            opts = optimoptions('fmincon', ...
                'Algorithm', 'interior-point', ...
                'InitBarrierParam', 100000, ...
                'ScaleProblem', true, ...
                'FiniteDifferenceType', 'central', ...
                'Display', 'off');
            
            if opt == 1
                algoparamsArray = repmat(struct('vinit', [], 'Sinit', [], 'sigEst', [], ...
                'lb', [], 'ub', [], 'solver', solver, 'options', opts), nvoxels, 1);

            elseif opt == 2
                algoparamsArray = repmat(struct('vinit', [], 'Sinit', [], 'sigEst', [], ...
                'lb', [], 'ub', [], 'solver', solver, 'options', opts), nvoxels, 1);

            elseif opt == 3
                 algoparamsArray = repmat(struct('vinit', [], 'Sinit', [], 'sigEst', [], ...
                'lb', [], 'ub', [], 'FixedWl', [], 'FixedR2wl', [], ...
                'solver', solver, 'options', opts), nvoxels, 1);

            elseif opt == 4
                algoparamsArray = repmat(struct('vinit', [], 'Sinit', [], 'sigEst', [], ...
                'lb', [], 'ub', [], 'solver', solver, 'options', opts), nvoxels, 1);

            elseif opt == 5
                algoparamsArray = repmat(struct('vinit', [], 'Sinit', [], 'sigEst', [], ...
                'lb', [], 'ub', [], 'FixedWl', [], 'FixedR2wl', [], ...
                'solver', solver, 'options', opts), nvoxels, 1);

            else
                error('Unknown option value: %d', opt);
            end

            
            vinit = 0.1;    % initial guess for R2 (for R2fat and R2water_long)
            %vmax = 5;       % maximum R2 value (ms^-1)
            vmax = 0.8;       % maximum R2 value (ms^-1)
            vmin = 0.5e-3;  % minimum R2 value (ms^-1)
            [maxS, idx] = max(S, [], 2);
            Sinit = maxS .* exp(vinit * echotimes(idx));
            
            if opt == 1
                % For 3-component model, split water equally initially.
                parfor i = 1:nvoxels
                    % Use half-and-half for water initialization.
                    algoparamsArray(i).vinit = vinit;
                    algoparamsArray(i).Sinit = Sinit(i);
                    algoparamsArray(i).sigEst = sigmaEstimateFromRoi;
                    algoparamsArray(i).lb = [0; 0; 0; 0.05; 0.249];
                    algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); 3 * Sinit(i);  0.25; 100];
                end
            elseif opt == 2
                % 2-component mode for sigma...
                parfor i = 1:nvoxels
                    algoparamsArray(i).vinit = vinit;
                    algoparamsArray(i).Sinit = Sinit(i);
                    algoparamsArray(i).sigEst = sigmaEstimateFromRoi;
                    algoparamsArray(i).lb = [0; 0; vmin; vmin; sigmaEstimateFromRoi*0.2];
                    algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); vmax; vmax; sigmaEstimateFromRoi*2];
                end
            elseif opt == 3
                % reduced model: [F; Ws; R2fat; R2water_short]
                if isempty(fixedWl_roi) || isempty(fixedR2wl_roi)
                    error('Fixed long-water maps required for opt==3');
                end
                parfor i = 1:nvoxels
                    algoparamsArray(i).vinit = vinit;
                    algoparamsArray(i).Sinit = Sinit(i);
                    algoparamsArray(i).sigEst = sigmaEstimateFromRoi;
                    % Lower / upper bounds for reduced fit:
                    % F and Ws cannot exceed 3*Sinit (rough)
                    % R2fat & R2water_short bounds keep same ranges

                    algoparamsArray(i).lb = [0;  vmin];  % enforce R2water_short >= 0.249 (T2* <= 4 ms) as before-ish
                    algoparamsArray(i).ub = [3* Sinit(i);  vmax];
                   
                    % algoparamsArray(i).lb = [0; 0.249];  % enforce R2water_short >= 0.249 (T2* <= 4 ms) as before-ish
                    % algoparamsArray(i).ub = [3*Sinit(i); 100];
                    algoparamsArray(i).FixedWl = fixedWl_roi(i);
                    algoparamsArray(i).FixedR2wl = fixedR2wl_roi(i);
                end
            elseif opt == 4
                % 2-component model.
                parfor i = 1:nvoxels
                    algoparamsArray(i).vinit = vinit;
                    algoparamsArray(i).Sinit = Sinit(i);
                    algoparamsArray(i).sigEst = sigmaEstimateFromRoi;
                    algoparamsArray(i).lb = [0; 0; vmin; vmin];
                    algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); vmax; vmax];
                end
            elseif opt == 5
                parfor i = 1:nvoxels
                    algoparamsArray(i).vinit = vinit;
                    algoparamsArray(i).Sinit = Sinit(i);
                    algoparamsArray(i).sigEst = sigmaEstimateFromRoi;
                    algoparamsArray(i).lb = [ ...
                    0;          % F
                    0;          % Ws
                    vmin;         % R2fat   (T2* ≤ 50 ms)
                    vmin];       % R2short (T2* ≤ 5 ms)

                    algoparamsArray(i).ub = [ ...
                    1.5*Sinit(i);
                    1.5*Sinit(i);
                    vmax;
                    100];

                    algoparamsArray(i).FixedWl   = fixedWl_roi(i);
                    algoparamsArray(i).FixedR2wl = fixedR2wl_roi(i);
                end
            
            else
                error('Unknown option value: %d. Valid options are 1,2,3 or 4', opt);
            end
    
        end 

        function results = RicianMagnitudeFitting(echotimes, tesla, Smagnitude, sig, algoparams)
            % unchanged full 6-parameter fitter (your original)
            RicianFit.solver = algoparams.solver;
            RicianFit.options = algoparams.options;
            RicianFit.lb = algoparams.lb;
            RicianFit.ub = algoparams.ub;
            RicianFit.objective = @(p) -FittingProcessorRician_3Comp.R2RicianObj(p, echotimes, tesla, Smagnitude, sig);
            
            initialGuesses = { [0.001; algoparams.Sinit/2; algoparams.Sinit/2; algoparams.vinit; algoparams.vinit; 0.2], ...
                               [algoparams.Sinit; 0.001; 0.001; algoparams.vinit; algoparams.vinit; 0.2]};
            nInit = numel(initialGuesses);
            pminCell = cell(nInit,1);
            fminArray = zeros(nInit,1);
            for k = 1:nInit
                RicianFit.x0 = initialGuesses{k};
                [p, fval] = fmincon(RicianFit);
                pminCell{k} = p;
                fminArray(k) = fval;
            end
            [results.fmin, idx] = min(fminArray);
            results.chosenmin = idx;
            bestP = pminCell{results.chosenmin};
            results.F = bestP(1);
            results.Wl = bestP(2);
            results.Ws = bestP(3);
            results.R2fat = bestP(4);
            results.R2water_long = bestP(5);
            results.R2water_short = bestP(6);
            results.SSE = FittingProcessorRician_3Comp.GetSEE(bestP, echotimes, tesla, Smagnitude);
        end


        %%%%%%%%% another reduced fit this time 3 exponeentila fit with
        %%%%%%%%% muscle para known%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function results = RicianMagnitudeFitting_FixMuscle( ...
            echotimes, tesla, Smagnitude, sig, algoparams)

          % Fixed (known)
         Wl   = algoparams.FixedWl;
         R2wl = algoparams.FixedR2wl;

        % Free parameters:
        % p = [F; Ws; R2fat; R2water_short]

        RicianFit.solver  = algoparams.solver;
        RicianFit.options = algoparams.options;
        RicianFit.lb = algoparams.lb;
        RicianFit.ub = algoparams.ub;

        RicianFit.objective = @(p) - ...
        FittingProcessorRician_3Comp.R2RicianObj_FixMuscle( ...
            p, echotimes, tesla, Smagnitude, sig, Wl, R2wl);

        % ---- Robust initializations ----
        S0 = algoparams.Sinit;

        p0_1 = [1.5*S0; 0.001; algoparams.vinit; algoparams.vinit];   % muscle-dominant
        p0_2 = [0.001;  1.5*S0; algoparams.vinit; algoparams.vinit]; % fat-dominant

        initialGuesses = {p0_1, p0_2};

        fmin = inf;
        bestP = [];

        for k = 1:numel(initialGuesses)
            RicianFit.x0 = initialGuesses{k};
            [p, fval] = fmincon(RicianFit);
                if fval < fmin
                fmin = fval;
                bestP = p;
                end
        end

        % ---- Output ----
        results.F  = bestP(1);
        results.Ws = bestP(2);
        results.R2fat = bestP(3);
        results.R2water_short = bestP(4);

        results.Wl = Wl;
        results.R2water_long = R2wl;
        results.S0 = results.F + results.Ws + Wl;
        results.fmin = fmin;
    end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% New reduced-fit function: fix long-water (Wl and R2wl) and fit remaining %%
        function results = RicianMagnitudeFitting_FixLong(echotimes, tesla, Smagnitude, sig, algoparams)
            % Parameterization:
            % Fixed: Wl = algoparams.FixedWl, R2wl = algoparams.FixedR2wl
            % Free vector: p = [F; Ws; R2fat; R2water_short]
            %
            % We optimize Rician loglik using fmincon with several initializations.
            
            % Create objective wrapper that expands p -> full parameter vector
            
            % Setup fmincon struct
            RicianFit.solver = algoparams.solver;
            RicianFit.options = algoparams.options;
            % Use reduced lb/ub from algoparams

            RicianFit.lb = algoparams.lb;
            RicianFit.ub = algoparams.ub;

            obj = @(p) -FittingProcessorRician_3Comp.R2RicianObj_FixLong(p, echotimes, tesla, Smagnitude, sig, algoparams.FixedWl, algoparams.FixedR2wl);
            RicianFit.objective = obj;

            % RicianFit.lb = [0;0];
            % RicianFit.ub = [3,200];
            
            % Initial guesses:

           % Ws_guess  = max( Smagnitude(1) - algoparams.FixedWl,  0.001 );
            %Ws_guess  = 0.001;
            %R2_guess  = algoparams.vinit;   % user-provided typical R2* guess       
            %R2_guess  = 1.5;   % user-provided typical R2* guess  
          %  R2_guess  = 15;   % user-provided typical R2* guess 
            % p0 guess 1: small fat, most signal assigned to Ws

            p0_1 = [algoparams.Sinit; algoparams.vinit];
            p0_2 = [0.001; algoparams.vinit];
           % p0_1 = [0.001; algoparams.Sinit - algoparams.FixedWl; algoparams.vinit; 0.3];
            % p0 guess 2: most signal assigned to fat (unlikely but try)
            %p0_2 = [algoparams.Sinit*0.5; algoparams.Sinit*0.25; algoparams.vinit; 0.3];
            initialGuesses = {p0_1, p0_2};
            
            nInit = numel(initialGuesses);
            pminCell = cell(nInit,1);
            fminArray = zeros(nInit,1);
            for k = 1:nInit
                RicianFit.x0 = initialGuesses{k};
                [p, fval] = fmincon(RicianFit);
                pminCell{k} = p;
                fminArray(k) = fval;
            end
            [fmin, idx] = min(fminArray);
            bestP = pminCell{idx};
            % Assign outputs
            results.fmin = fmin;
            results.chosenmin = idx;
            results.Wl = algoparams.FixedWl;
            results.Ws = bestP(1);
            results.R2water_long = algoparams.FixedR2wl;
            results.R2water_short = bestP(2);
            %results.SSE = FittingProcessorRician_3Comp.GetSEE([bestP(1); algoparams.FixedWl; bestP(2); bestP(3); algoparams.FixedR2wl; bestP(4)], echotimes, tesla, Smagnitude);
        end

%%%%RAndy%%%%%%%%%%%%%%%%%%%%
        function results = RicianMagnitudeFitting_WithSigma_randy(echotimes, tesla, Smagnitude, sig, algoparams)
            % unchanged
            Ricianfitting.solver = algoparams.solver;
            Ricianfitting.options = algoparams.options;

            Ricianfitting.lb = algoparams.lb;
            Ricianfitting.ub = algoparams.ub;


            Ricianfitting.objective = @(p) -FittingProcessorRician_3Comp.R2RicianObj_WithSigmarandy(p, echotimes, tesla, Smagnitude,sig);
            Ricianfitting.x0 = [0.001; algoparams.Sinit; algoparams.vinit; algoparams.vinit];
            [pmin1, fmin1] = fmincon(Ricianfitting);
            Ricianfitting.x0 = [algoparams.Sinit; 0.001; algoparams.vinit; algoparams.vinit];
            [pmin2, fmin2] = fmincon(Ricianfitting);
            
            if fmin1 <= fmin2
                results.F = pmin1(1);
                results.W = pmin1(2);
                results.R2water = pmin1(3);
                results.R2fat = pmin1(4);
                %results.sig = pmin1(5);
                results.fmin = fmin1;
                results.chosenmin = 1;
              %  results.SSE = FittingProcessorRician_3Comp.GetSEE(pmin1, echotimes, tesla, Smagnitude);
            else
                results.F = pmin2(1);
                results.W = pmin2(2);
                results.R2water = pmin2(3);
                results.R2fat = pmin2(4);
                %results.sig = pmin2(5);
                results.fmin = fmin2;
                results.chosenmin = 2;
               % results.SSE = FittingProcessorRician_3Comp.GetSEE(pmin2, echotimes, tesla, Smagnitude);
            end
          end


        %%%%% end of randy%%%%%%%%%%%%%%%%%%%%%%%%
        
        function results = RicianMagnitudeFitting_WithSigma(echotimes, tesla, Smagnitude, sigmaInit, algoparams)
            % unchanged
            Ricianfitting.solver = algoparams.solver;
            Ricianfitting.options = algoparams.options;
            Ricianfitting.lb = algoparams.lb;
            Ricianfitting.ub = algoparams.ub;
            Ricianfitting.objective = @(p) -FittingProcessorRician_3Comp.R2RicianObj_WithSigma(p, echotimes, tesla, Smagnitude);

            Ricianfitting.x0 = [0.001; algoparams.Sinit; algoparams.vinit; algoparams.vinit; sigmaInit];
            [pmin1, fmin1] = fmincon(Ricianfitting);
            Ricianfitting.x0 = [algoparams.Sinit; 0.001; algoparams.vinit; algoparams.vinit; sigmaInit];
            [pmin2, fmin2] = fmincon(Ricianfitting);

            if fmin1 <= fmin2
                results.F = pmin1(1);
                results.W = pmin1(2);
                results.R2water = pmin1(3);
                results.R2fat = pmin1(4);
                results.sig = pmin1(5);
                results.fmin = fmin1;
                results.chosenmin = 1;
                results.SSE = FittingProcessorRician_3Comp.GetSEE(pmin1, echotimes, tesla, Smagnitude);
            else
                results.F = pmin2(1);
                results.W = pmin2(2);
                results.R2water = pmin2(3);
                results.R2fat = pmin2(4);
                results.sig = pmin2(5);
                results.fmin = fmin2;
                results.chosenmin = 2;
                results.SSE = FittingProcessorRician_3Comp.GetSEE(pmin2, echotimes, tesla, Smagnitude);
            end
        end
        %%%%%%%%%%%%%%%%expreiment done %%%%%%%%%%%%%

        function sse = GetSEE(p, echotimes, tesla, Smeasured)
            % unchanged
            if numel(p) == 6
                Spredicted = FittingProcessorRician_3Comp.MultiPeakFatTripleR2(echotimes, tesla, ...
                    p(1), p(2), p(3), p(4), p(5), p(6), 0);
            elseif numel(p) == 5
                Spredicted = FittingProcessorRician_3Comp.MultiPeakFatDoubleR2(echotimes, tesla, ...
                    p(1), p(2), p(3), p(4), 0);
            else
                error('Parameter vector has insufficient elements.');
            end
            errors = abs(Spredicted) - abs(Smeasured);
            sse = sum(errors.^2);
        end

        
        function [loglik] = R2RicianObj(p, echotimes, tesla, Smeasured, sig)
            % unchanged
            Spredicted = abs(FittingProcessorRician_3Comp.MultiPeakFatTripleR2(echotimes, tesla, p(1), p(2), p(3), p(4), p(5), p(6),0));
            Smeasured = abs(Smeasured);
            loglik = FittingProcessorRician_3Comp.RicianLogLik(Smeasured, Spredicted, sig);
        end

        %%% New wrapper for reduced-fit objective %%%
        function [loglik] = R2RicianObj_FixLong(p, echotimes, ~, Smeasured, sig, FixedWl, FixedR2wl)%tesla
            % p = [F; Ws; R2fat; R2water_short]
            % Expand to full vector expected by MultiPeakFatTripleR2:
            % [F; Wl; Ws; R2fat; R2water_long; R2water_short]
           % fullP = [p(1); FixedWl; p(2); p(3); FixedR2wl; p(4)];
           %F_short = p(1);
           %R2_short = p(2);
            %fullP = [p(1); FixedWl; p(2); FixedR2wl];
            %Spredicted = abs(FittingProcessorRician_3Comp.MultiPeakFatTripleR2(echotimes, tesla, fullP(1), fullP(2), fullP(3), fullP(4), fullP(5), fullP(6), 0));

            Spredicted = abs(FittingProcessorRician_3Comp.TwoWaterDoubleR2(echotimes, p(1), FixedWl,FixedR2wl, p(2), 0));
            %Spredicted = abs(FittingProcessorRician_3Comp.MultiPeakFatDoubleR2(echotimes, tesla, p(1), FixedWl, FixedR2wl, p(2), 0));
            Smeasured = abs(Smeasured);
            loglik = FittingProcessorRician_3Comp.RicianLogLik(Smeasured, Spredicted, sig);
        end

        %%%%%%%%%%%%%%%Fix Muscle objective%%%%%%%%%%%%%%%%%

        function loglik = R2RicianObj_FixMuscle( ...
            p, echotimes, tesla, Smeasured, sig, Wl, R2wl)

            % p = [F; Ws; R2fat; R2water_short]

                F   = p(1);
                Ws  = p(2);
                R2f = p(3);
                R2s = p(4);

                Spred = abs( ...
                    FittingProcessorRician_3Comp.MultiPeakFatTripleR2( ...
                    echotimes, tesla, ...
                    F, Wl, Ws, ...
                    R2f, R2wl, R2s, 0) );

                loglik = FittingProcessorRician_3Comp.RicianLogLik( ...
                abs(Smeasured), Spred, sig);
        end

        %%%%%%%%%%%%%%%Fix Muscle objective end %%%%%%%%%%%%%%%%%%%%%%%

        function [loglik] = R2RicianObj_WithSigma(p, echotimes, ~, Smeasured)
            % unchanged
            %Spredicted = abs(FittingProcessorRician_3Comp.MultiPeakFatDoubleR2(echotimes, tesla, p(1), p(2), p(3), p(4), 0));
            Spredicted = abs(FittingProcessorRician_3Comp.TwoWaterDoubleR2(echotimes, p(1), p(2), p(3), p(4), 0));
            Smeasured = abs(Smeasured);
            loglik = FittingProcessorRician_3Comp.RicianLogLik(Smeasured, Spredicted, p(5));
        end

        function [loglik] = R2RicianObj_WithSigmarandy(p, echotimes, tesla, Smeasured, sig)
            % unchanged
            Spredicted = abs(FittingProcessorRician_3Comp.MultiPeakFatDoubleR2(echotimes, tesla, p(1), p(2), p(3), p(4), 0));
            %Spredicted = abs(FittingProcessorRician_3Comp.TwoWaterDoubleR2(echotimes, p(1), p(2), p(3), p(4), 0));
            Smeasured = abs(Smeasured);
            loglik = FittingProcessorRician_3Comp.RicianLogLik(Smeasured, Spredicted, sig);
        end
            
        
        function S = MultiPeakFatTripleR2(t, tesla, F, Wl, Ws, R2fat, R2water_long, R2water_short, fB)
            % unchanged
            gyro = 42.58e6;
            larmor = tesla * gyro;
            %Fatamps = [0.087; 0.694; 0.128; 0.004; 0.039; 0.048];
            Fatamps  =[0.089; 0.577; 0.059; 0.093; 0.059; 0.013; 0.02; 0.02; 0.01; 0.059];
            Fatshift =[-3.81; -3.4; -3.12; -2.67; -2.45; -1.94; -0.63; -0.4; 0.52; 0.62];
            %Fatshift = [-3.9; -3.5; -2.7; -2.04; -0.49; 0.50];
            FatFreqs = Fatshift * larmor * 1e-6;
            FatFreqs = FatFreqs / 1000;
            Fatw = 2*pi*FatFreqs;
            Waterw = 0;
            fatExp = exp(1i * Fatw * t');  % 6-by-T
            fatSignal = Fatamps.' * fatExp;  % 1-by-T
            waterExp = exp(1i * Waterw * t');  % 1-by-T
            fatComponent = F * fatSignal .* exp(-R2fat*t');
            waterLongComponent = Wl * waterExp .* exp(-R2water_long*t');
            waterShortComponent = Ws * waterExp .* exp(-R2water_short*t');
            S = fatComponent + waterLongComponent + waterShortComponent;
            fB = fB/1000;
            offResonance = exp(1i * 2*pi*fB*t');
            S = S .* offResonance;
        end

        % function S = MultiPeakFatTripleR2(t, tesla, F, Wl, Ws, ...
        %                          R2water_long, R2water_short, fB)
        % 
        %     gyro   = 42.58e6;
        %     larmor = tesla * gyro;
        % 
        %     % ----- Fat model -----
        %     Fatamps  = [0.089; 0.577; 0.059; 0.093; 0.059; ...
        %                 0.013; 0.02; 0.02; 0.01; 0.059];
        % 
        %     Fatshift = [-3.81; -3.4; -3.12; -2.67; -2.45; ...
        %                  -1.94; -0.63; -0.4; 0.52; 0.62];
        % 
        %     FatR2 = [130 135 85 95 100 90 100 100 85 85];   % s^-1
        %     FatR2 = FatR2(:);                               % column vector
        % 
        %     % Frequencies (Hz → kHz → rad/s)
        %     FatFreqs = Fatshift * larmor * 1e-6 / 1000;
        %     Fatw     = 2*pi*FatFreqs;
        % 
        %     % Time vector
        %     t = t(:).';   % 1 × T
        % 
        %     % Fat signal with per-peak R2*
        %     fatExp = exp(1i * Fatw * t) .* exp(-FatR2 * t);   % (10 × T)
        %     fatSignal = Fatamps.' * fatExp;                   % (1 × T)
        %     fatComponent = F * fatSignal;
        % 
        %     % ----- Water components -----
        %     waterExp = ones(size(t));  % exp(i*0*t)
        % 
        %     waterLongComponent  = Wl * waterExp .* exp(-R2water_long  * t);
        %     waterShortComponent = Ws * waterExp .* exp(-R2water_short * t);
        % 
        %      % ----- Combine -----
        %     S = fatComponent + waterLongComponent + waterShortComponent;
        % 
        %     % ----- Off-resonance -----
        %     fB = fB / 1000;  % Hz → kHz
        %     offResonance = exp(1i * 2*pi*fB*t);
        %     S = S .* offResonance;
        % end

        
        function S = MultiPeakFatDoubleR2(t, tesla, F, W, R2water, R2fat, fB)
            % unchanged
            gyro = 42.58e6;
            larmor = tesla * gyro;

            % old 6 point

            % Fatamps = [0.087; 0.694; 0.128; 0.004; 0.039; 0.048];
            % Fatshift = [-3.9; -3.5; -2.7; -2.04; -0.49; 0.50];

            % new 10 point one 
            Fatamps = [0.089; 0.577; 0.059; 0.093; 0.059; 0.013; 0.02; 0.02; 0.01; 0.059];
            Fatshift = [-3.81; -3.4; -3.12; -2.67; -2.45; -1.94; -0.63; -0.4; 0.52; 0.62];

            FatFreqs = Fatshift * larmor * 1e-6;
            FatFreqs = FatFreqs/1000;
            Fatw = 2*pi*FatFreqs;
            Waterw = 0;
            fatExp = exp(1i * Fatw * t');
            fatSignal = Fatamps.' * fatExp;
            waterExp = exp(1i * Waterw * t');
            fatComponent = F * fatSignal .* exp(-R2fat*t');
            waterComponent = W * waterExp .* exp(-R2water*t');
            S = fatComponent + waterComponent;
            fB = fB/1000;
            offResonance = exp(1i * 2*pi*fB*t');
            S = S .* offResonance;
        end

        function S = TwoWaterDoubleR2(t, F, W, R2water, R2fat, fB)

             % t must be column or row (ms or seconds — keep consistent)
             t = t(:).';                     % force row vector

             % nEchoUse = min(4, numel(t));
             % t = t(1:nEchoUse);

            % Known long-T2* component
            longComp  = W .* exp(-R2water .* t);

             % Unknown short-T2* component
            shortComp = F .* exp(-R2fat .* t);

             % Combine
             S = shortComp + longComp;

            % Off-resonance term (fB assumed in Hz)
             freqs = exp(1i * 2*pi*(fB/1000) .* t);   % divide by 1000 if your old model used kHz
              S = S .* freqs;

        end

        % function S = TwoWaterDoubleR2_new(t, F, W, R2water, R2fat, fB)
        % 
        %      % t must be column or row (ms or seconds — keep consistent)
        %      t = t(:).';                     % force row vector
        % 
        %       nEchoUse = min(4, numel(t));
        %       t = t(1:nEchoUse);
        % 
        %     % Known long-T2* component
        %     longComp  = W .* exp(-R2water .* t);
        % 
        %      % Unknown short-T2* component
        %     shortComp = F .* exp(-R2fat .* t);
        % 
        %      % Combine
        %      S = shortComp + longComp;
        % 
        %     % Off-resonance term (fB assumed in Hz)
        %      freqs = exp(1i * 2*pi*(fB/1000) .* t);   % divide by 1000 if your old model used kHz
        %       S = S .* freqs;
        % 
        % end


        
        function [loglik, logliks] = RicianLogLik(measurements, predictions, sigma)
            sigmaSquared = sigma.^2;
            sumsqsc = (measurements.^2 + predictions.^2) ./ (2*sigmaSquared);
            scp = measurements.*predictions ./ sigmaSquared;
            lb0 = FittingProcessorRician_3Comp.logbesseli0(scp);
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
        
        function loglik = GaussianLogLik(measured, predicted, sigma)
            errors = measured - predicted;
            n = numel(errors);
            loglik = -0.5*sum((errors./sigma).^2) - n*log(sigma*sqrt(2*pi));
        end
        
        function img_whiten = WhitenImage(img)
            img_whiten = img;
            num_slices = size(img,3);
            for s = 1:num_slices
                slice = img(:,:,s);
                centered = slice - median(slice(:));
                normalized = centered / prctile(centered(:),75);
                p25 = prctile(normalized(:),25);
                p75 = prctile(normalized(:),75);
                scaleFactor = 2 / (p75 - p25);
                offset = -1 - 2*p25/(p75-p25);
                img_clean = normalized*scaleFactor + offset;
                lower_bound = prctile(img_clean(:),3);
                upper_bound = prctile(img_clean(:),97);
                img_clean = max(min(img_clean,upper_bound),lower_bound);
                img_whiten(:,:,s) = img_clean;
            end
        end
    end
end
