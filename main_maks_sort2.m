clc; clear; close all;

%% ===== Load Multi-label NIfTI =====
labelFile = 'water_reg_dseg.nii.gz';
labelImg  = niftiread(labelFile);
info      = niftiinfo(labelFile);
info.Datatype = 'uint8';

%% ===== Label → Muscle Name Mapping =====
labelMap = containers.Map( ...
    [7102, 7101, 7112, 7111, 7122, 7121, 7132, ...
     7131, 7162, 7161, 7172, 7171, 7182, 7181, 7192, 7191], ...
    { ...
     'VL-R', ...
     'VL-L', ...
     'VI-R', ...
     'VI-L', ...
     'VM-R', ...
     'VM-L', ...
     'RF-R', ...
     'RF-L', ...
     'SM-R', ...
     'SM-L', ...
     'ST-R', ...
     'ST-L', ...
     'BF-R', ...
     'BF-L', ...
     'BFS-R', ...
     'BFS-L' ...
    } ...
);

%% ===== Erosion Parameters (Small) =====
se    = strel('disk',1);   % gentle erosion
nIter = 2;

%% ===== Initialize Regional Masks =====
leftQuadMask  = false(size(labelImg));
rightQuadMask = false(size(labelImg));
leftHamsMask  = false(size(labelImg));
rightHamsMask = false(size(labelImg));

%% ===== Process Each Muscle =====
keys = cell2mat(labelMap.keys);

for i = 1:numel(keys)

    lbl = keys(i);

    if ~any(labelImg(:) == lbl)
        fprintf('Label %d not found → skipped\n', lbl);
        continue;
    end

    % ---- Binary mask ----
    mask = (labelImg == lbl);

    % ---- Slice-wise 2D erosion ----
    maskEroded = false(size(mask));
    for z = 1:size(mask,3)
        temp = mask(:,:,z);
        for k = 1:nIter
            temp = imerode(temp, se);
        end
        maskEroded(:,:,z) = temp;
    end

    % ---- Save individual muscle ----
    outName = sprintf('%s.nii', labelMap(lbl));
    niftiwrite(uint8(maskEroded), outName, info);
    fprintf('Saved %s\n', outName);

    %% ===== Add to Regional Masks =====
    switch labelMap(lbl)
        % ---- Quadriceps ----
        case {'VL-L','VI-L','VM-L','RF-L'}
            leftQuadMask  = leftQuadMask  | maskEroded;

        case {'VL-R','VI-R','VM-R','RF-R'}
            rightQuadMask = rightQuadMask | maskEroded;

        % ---- Hamstrings ----
        case {'SM-L','ST-L','BF-L','BFS-L'}
            leftHamsMask  = leftHamsMask  | maskEroded;

        case {'SM-R','ST-R','BF-R','BFS-R'}
            rightHamsMask = rightHamsMask | maskEroded;
    end
end

%% ===== Save Regional Masks =====
niftiwrite(uint8(leftQuadMask),  'leftQuadMask.nii',  info);
niftiwrite(uint8(rightQuadMask), 'rightQuadMask.nii', info);
niftiwrite(uint8(leftHamsMask),  'leftHamsMask.nii',  info);
niftiwrite(uint8(rightHamsMask), 'rightHamsMask.nii', info);

disp('===================================================');
disp('All muscle and regional masks eroded and saved');
disp('===================================================');
