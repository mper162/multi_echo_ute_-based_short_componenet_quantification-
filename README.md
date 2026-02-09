UTE Muscle Fat and Short/Long T2* Quantification Pipeline

This repository contains a multi-echo UTE MRI processing pipeline for quantifying muscle fat fraction, long T2* water, and short T2* components at individual muscle and regional levels.

The workflow integrates Dixon water–fat separation, fat subtraction, mono- and bi-exponential fitting, muscle segmentation, and batch Excel reporting.

UTE magnitude + phase
        ↓
Water–fat separation (late echoes)
        ↓
Fat subtraction from original UTE
        ↓
Mono-exponential fit (long T2*)
        ↓
Bi-exponential Rician fit (short + long T2*)
        ↓
Mask-based quantification
        ↓
Batch Excel output


Step-by-Step Workflow
1. Multi-echo UTE data preparation

Magnitude and phase UTE images are merged.

Input data consist of multi-echo UTE acquisitions.

2. Water–fat separation (Dixon-based)

Water–fat separation is performed using later echoes only.

This ensures short T2* components do not contaminate the fat estimation.

Script:

script_dixon_ute_siemens_final2 (Mathematica)

but source and output in a seperate space. this might be a problem extracting masks. so use 
fslcpgeom merged_ute_mag.nii.gz UTE_mag_original_0.nii after the ute_folder to convert generated outputs to the same space as input _UTE. Do it for every file.

3. Fat signal removal

The estimated fat signal is removed from the original UTE dataset.

Produces a fat-subtracted UTE image.

Script:

test_fat_sub_2

4. Mono-exponential fitting (long T2*)

Mono-exponential decay model applied to fat-subtracted data.

Uses later echoes only.

Outputs:

Long-component S₀

Long-component T2*

Purpose:

Robust estimation of muscle long-T2* signal.

Script:

Mono-exponential fitting module

5. Bi-exponential fitting (short + long T2*)

Bi-exponential model separating:

Short T2* component

Long T2* component

Rician noise model is used.

Outputs:

Short-component S₀

Relative short vs. long signal fractions

Script:

rician_randika3

6. Muscle fat and water fraction calculation (single mask)

Computes fat, long-water, and short-water fractions within a mask.

Script:

hist_total_muscle_fat_fraction

Output metrics:

Fat percentage

Long-water percentage

Short-water percentage

Final combined fractions

7. Muscle mask generation and labeling

Muscle masks are derived from Dixon muscle scans.

Workflow:

Dixon muscle segmentation

Registration between Dixon and UTE

Transfer of muscle labels to UTE space

Script:

main_mask_sort2

Outputs:

Individual muscle masks

Regional masks (quadriceps, hamstrings, left/right)

8. Batch processing and Excel reporting

All individual and regional masks are processed automatically.

For each mask, the following are computed:

Fat percentage

Total water percentage

Long-water percentage

Short-water percentage

Final combined fractions

Output:

Single Excel file summarizing all muscles and regions

Script:

final_excel_mask_percentages

9. Sanity check using synthetic 4D signal decay

T2* and S₀ maps are used to synthesize 4D signals across echo times.

Separate synthetic signals are generated for:

Short-T2* component

Long-T2* component

Expected behavior:

Short-T2* signal decays within the first few echoes

In practice, short components vanish after ~9 echoes (~2–3 ms), which is physiologically plausible

Script:

shortlong_4d_maps

Outputs

Muscle-specific fat, long-water, and short-water fractions

Regional quadriceps and hamstrings metrics

Excel summary table for statistical analysis

Optional synthetic 4D signal data for validation
