# Output files and manuscript mapping

Files below are created from the supplied input. No numerical result is bundled
in this repository. The run-specific `audit/output_inventory.json` records every
generated file, its SHA-256, byte count, and exact CSV column names.

## Core tables and plotting data

| File | Row unit | Main columns beyond identifying keys |
|---|---|---|
| `tables/unit_manifest.csv` | unit-session | Group, Mouse_ID, Session_ID, Region, Unit_ID, Stage |
| `tables/event_bins.csv` | recorded trial-bin shared by units | Group, Mouse_ID, Session_ID, Stage, Trial_ID, Trial_Order, Time_Bin_Start_s, Time_Bin_End_s, Sound, Light, Shock |
| `tables/unit_count_qc.csv` | unit-session | N_Trials, N_Bins, Total_Spikes, Exposure_Total_s, Zero_Bin_Fraction, Max_Count, Mean_FR_Hz |
| `tables/sample_allocation.csv` | mouse-session-region | N_Unit_Sessions, Min_Trials, Max_Trials, N_Bins, Total_Spikes |
| `tables/baseline_parameters.csv` | unit-session | Baseline_FR_Hz, Baseline_SD_Hz |
| `plotdata/unit_psth.csv` | unit-session-bin | Time_Bin_Start_s, Time_Bin_End_s, FR_Hz, N_Trials, Baseline_FR_Hz, Baseline_SD_Hz, Response_Z |
| `plotdata/animal_psth.csv` | mouse-region-stage-bin | FR_Hz, Response_Z, N_Units, N_Z_Units |
| `plotdata/group_psth.csv` | group-region-stage-bin | Mean_FR_Hz, SEM_FR_Hz, Mean_Z, SEM_Z, N_Animals, N_Z_Animals |
| `plotdata/unit_window_rates.csv` | unit-session-window | Event_FR_Hz, Baseline_FR_Hz, Delta_FR_Hz, N_Trials, Window, Window_Start_s, Window_End_s |
| `plotdata/animal_window_rates.csv` | mouse-region-stage-window | Event_FR_Hz, Baseline_FR_Hz, Delta_FR_Hz, N_Units, Min_Trials, Max_Trials, Log_Rate_Ratio |
| `statistics/response_summaries.csv` | group-region-stage-window-metric | Metric, N_Animals, Mean, CI_Lower, CI_Upper |
| `statistics/response_tests.csv` | planned response hypothesis | Family, Region, Window, Metric, P_Raw, P_Adjusted, Correction, Family_Size, Test, Assignments, effect and CI columns |

Unit identifiers are always scoped to session and region. Animal keys include
Group, Mouse_ID, Region and, where applicable, Stage. Unit-session tables are not
independent animal tables.

Paired statistics include N_Animals, Mean, CI_Lower, CI_Upper, Effect_dz. Group
comparisons include Group_A, Group_B, N_A, N_B, Mean_A, Mean_B,
Difference_A_Minus_B, CI_Lower, CI_Upper and Hedges_g. Means and CIs remain
unadjusted for multiplicity; only the probability columns receive correction.

## GLM outputs

| File | Row unit | Main columns |
|---|---|---|
| `tables/glm_trial_splits.csv` | unit-session-trial-fold | unit keys, Fold, Trial_ID, Trial_Order, Role |
| `tables/glm_folds.csv` | unit-session-fold | N_Train_Trials, N_Test_Trials, Heldout_Spikes, LL_Full_Constant_Omitted, LL_Reduced_Constant_Omitted, Delta_LL_Nats, Gain_Bits_Per_Spike |
| `tables/glm_coefficients.csv` | unit-session-fold-model-term | Model, Term, Beta, Exp_Beta, Training_Order_Mean, Training_Order_SD |
| `tables/glm_diagnostics.csv` | unit-session-fold | Full_Iterations, Reduced_Iterations, Train_Pearson_Dispersion, N_Train_Bins, Status |
| `plotdata/glm_unit_gains.csv` | unit-session-term | Term, Gain_Bits_Per_Spike, Delta_LL_Nats, Heldout_Spikes, N_Trials |
| `plotdata/glm_animal_gains.csv` | mouse-region-stage-term | Term, Gain_Bits_Per_Spike, N_Units |
| `plotdata/glm_unit_predictions.csv` | unit-session-bin | Time_Bin_Start_s, Observed_FR_Hz, Predicted_Full_FR_Hz, Predicted_Reduced_FR_Hz |
| `plotdata/glm_animal_predictions.csv` | mouse-region-stage-bin | same prediction columns |
| `plotdata/glm_animal_null.csv` | mouse-region-stage-draw | Permutation_ID, Gain_Bits_Per_Spike |
| `plotdata/glm_group_null.csv` | group-region-stage-draw | Permutation_ID, Gain_Bits_Per_Spike |
| `statistics/glm_event_locking.csv` | planned global event test | Group, Region, Stage, Term, N_Animals, N_Units, Mean, CI_Lower, CI_Upper, N_Shifts, P_Raw, P_Adjusted, Correction, Family_Size, Family, Test |
| `statistics/glm_drop_one_descriptive.csv` | group-region-stage-term | N_Animals, Mean, CI_Lower, CI_Upper |

## RSA outputs

| File | Row unit | Main columns |
|---|---|---|
| `plotdata/rsa_profiles.csv` | mouse-region-condition-split-feature | Group, Mouse_ID, Region, Stage, Feature, Time_s, Value, N_Valid_Units, N_Total_Units, Min_Trials, Max_Trials, Condition, N_Features_Expected, Split |
| `plotdata/rsa_matrices.csv` | mouse-region-directed matrix edge | Condition_1, Condition_2, Similarity_r, N_Features, Status |
| `plotdata/rsa_reliability.csv` | mouse-region-condition | Reliability_r, Odd_Min_Trials, Odd_Max_Trials, Odd_Min_Valid_Units, Even_Min_Trials, Even_Max_Trials, Even_Min_Valid_Units |
| `statistics/rsa_primary.csv` | planned region contrast | group-comparison columns, Family, Region, Condition_1, Condition_2, Mean_r_A, Mean_r_B, Test_Scale, N_Clipped_A, N_Clipped_B, P_Adjusted, Correction, Family_Size |
| `statistics/rsa_edge_summaries.csv` | group-region-edge | N_Animals, Mean, CI_Lower, CI_Upper on r scale |
| `plotdata/rsa_sensitivity.csv` | mouse-region-variant | directed primary edge columns, Variant |
| `statistics/rsa_sensitivity.csv` | region-variant hypothesis, Fig. 2 only | primary group-test columns, Variant, P_Adjusted, Correction, Family_Size |
| `plotdata/rsa_count_control_repeats.csv` | mouse-region-variant-repetition | Repeat, Keep_Probability, N_Selected_Units, N_CL_Valid_Units, N_TS_Valid_Units, Similarity_r |
| `plotdata/rsa_count_controls.csv` | mouse-region-variant | Keep_Probability, Median_r, Monte_Carlo_Lower, Monte_Carlo_Upper, N_Repeats, N_Finite |
| `statistics/rsa_count_controls_descriptive.csv` | group-region-variant | N_Animals, Mean, CI_Lower, CI_Upper of animal medians |

Diagonal RSA entries are self-correlations of estimable profiles, not split-half
reliability. Missing profiles remain undefined. Temporal-profile correlations
do not compare the identity of individual cells.

## Optional behavior

`plotdata/behavior_animals.csv` retains Mouse_ID, Group, Day, Freezing_Percent.
`plotdata/behavior_groups.csv` contains Group, Day, Mean, SEM, N_Animals.
`statistics/behavior_tests.csv` contains Day, Group_A, Group_B, N_A, N_B,
Difference_A_Minus_B, T, DF, P_Raw, P_Adjusted, Correction, Family_Size, Family, Test.

## Figure mapping

Main panels are exported as `Fig2` (with optional behavior) or
`Fig2_Electrophysiology_BCD`, `Fig4_Electrophysiology_JKL`, and `V1_Main`.
Each is exported as PDF, SVG, and PNG. Overall titles are omitted; panel letters,
region/condition labels, axis units, and color scales are retained.

| Figure output stem | Evidence shown | Source |
|---|---|---|
| Fig. 2B | Inference heatmaps and animal-weighted z PSTHs | unit_psth, group_psth |
| Fig. 2C | Global event predictive gains | glm_animal_gains, glm_event_locking |
| Fig. 2D | Shared CL–TS matrices and animal group contrast | rsa_matrices, rsa_primary |
| Fig. 4J | Sound rate by group and region | animal_window_rates, response_tests |
| Fig. 4K | Sound/baseline log ratio | same |
| Fig. 4L | Three-condition matrices and CL–TS contrast | rsa_matrices, rsa_primary |
| V1 A–C | Unit heatmap, rate PSTH, paired sound response | unit_psth, group_psth, animal_window_rates, response_tests |
| Supplement_Control_Population | No pair available-stage heatmaps | unit_psth, group_psth |
| Supplement_Count_QC | Units/mouse, zero bins, baseline rates | sample_allocation, unit_count_qc, baseline_parameters |
| Supplement_Raw_PSTH | Unstandardized group time courses | group_psth |
| Supplement_Response_Windows | Animal response increments | animal_window_rates |
| Supplement_Intervention_Rates | Baseline rate and absolute change | animal_window_rates, response_tests |
| Supplement_GLM_DropOne | Descriptive event contributions | glm_animal_gains |
| Supplement_GLM_Predictions | Held-out observed and predicted rates | glm_animal_predictions |
| Supplement_GLM_Nulls | Shift-null distributions and observed statistic | glm_group_null, glm_event_locking |
| Supplement_GLM_Diagnostics | Dispersion and convergence | glm_diagnostics |
| Supplement_RSA_Profiles | Actual time courses used as features | rsa_profiles |
| Supplement_RSA_Available_Matrices | Additional Inference three-condition matrix | rsa_matrices |
| Supplement_RSA_Reliability | Odd/even profile correlations | rsa_reliability |
| Supplement_RSA_Sensitivity | All configured definitions | rsa_sensitivity |
| Supplement_RSA_CountControls | Animal medians over count controls | rsa_count_controls |

Supplementary output names are semantic rather than hard-coded global S numbers,
so they can be assigned manuscript S numbers after integrating other experiments.
Fig. 2 and Fig. 4 directories remain separate to prevent filename collisions.
Tables use semantic names for the same reason; their data are ready for final
supplementary table assembly.

## Audit

- `run_manifest.json`: completion state, code/package version, Git commit when
  available, Python/dependency versions, configuration and source hashes.
- `input_summary.json`: actual animal/session/unit-session/trial/bin counts.
- `exclusions.csv`: deliberate stage omissions and invalid unit-trial epochs.
- `normalization_unavailable.csv`, `glm_failures.csv`, `rsa_unavailable.csv`,
  `count_control_exclusions.csv`: analysis-specific availability.
- `config_used.json`, `figure_settings.json`: exact settings.
- `output_inventory.json`: actual output filenames, hashes and column schemas.

These audit outputs are part of reproducibility, not extra independent samples.
