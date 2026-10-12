# Figure legend information

Use final run outputs to supply all counts and probabilities. These are legend
templates, not prewritten results. Insert counts from `sample_allocation.csv`,
`input_summary.json`, and the analysis-specific tables; the number of analyzable
animals can differ by metric. Do not use a generic n for all sampling levels.

## Fig. 2

**A**, if supplied: session-level freezing, mean ± SEM across independent mice.
State animals/group, days, and the Welch/Holm family used. **B**, unit-session
baseline-z heatmaps arranged by animal and unit identifier, with animal boundaries;
rows are not matched across stages. Curves show animal means ± SEM. State
baseline window, bin width, actual stimulus intervals, and clipping limits for
color display. Heatmap clipping does not alter analysis values. Report animals,
recording sessions and unit-sessions per stage/region. **C**, animal mean held-out
GLM gain, with animal points and black mean/percentile CI; global circular-shift
tests and the BH correction family are reported in the statistical table. Positive
gain, a CI excluding zero, and shift significance are separate quantities.
**D**, group means of animal-level CL–TS matrices and animal correlations;
two-sided group-label tests on Fisher z with Holm correction across regions.
No pair pairing-stage recordings are outside the selected analysis. Define CL
and TS. Report correlation estimation windows and N used in each regional test.

## Fig. 4 electrophysiology

**J**, animal mean test-sound rate. **K**, natural logarithm of the animal
sound/pre-cue rate ratio, calculated after unit averaging. **L**, available
three-condition temporal-profile matrices and CL–TS animal correlations.
State independent group sizes, pretest treatment timing, and that both baseline
and sound were recorded after treatment. J/K probabilities belong to a shared
four-test family; regional RSA tests form a separate two-test Holm family.
Markers show mice and black diamonds show means with percentile animal CIs.
Reduced profile reliability is evaluated separately and can limit RSA interpretation.

## V1

**A**, sound-aligned unit-session heatmap with animal boundaries. **B**, raw-rate
PSTH averaged first across trials and units within mice, then across mice;
shading is animal SEM. **C**, paired animal baseline and full-sound rates,
with exact or Monte Carlo sign-flip probability as labeled by the actual run.
State animals, sessions, units, retained trials and windows. Early/late windows
are a separate two-test Holm family; light activity is descriptive. This design
has no pretraining or unpaired V1 comparison.

## Supplementary figures

- Count QC: small points represent units; emphasized summaries are animal-level.
  These are count characteristics, not waveform/ISI isolation metrics.
- GLM predictions: solid observed, dashed full model, dotted reduced model;
  colors identify groups. Curves contain only held-out observations/predictions.
- GLM nulls: bars are resampling draws, colored line is the observed animal-mean
  gain. Draws are not biological replicates. State shared shifts and draw count.
- GLM diagnostics: each observation is a unit-fold fit; these plots are descriptive.
- RSA reliability: odd/even baseline normalization is separate, and actual
  split trial counts appear in the source table.
- RSA sensitivity: state all windows, metrics and bin widths, including that
  500-ms profiles average already standardized 250-ms bins. Report the configured
  family; intervention variants are descriptive.
- Count controls: each point is an animal's median over repetitions. State the
  selected unit number, repetitions, retained animal count, and data-derived
  thinning probability. Monte Carlo repetition intervals and animal bootstrap
  intervals must not be confused. The procedure does not match reliability.

Mouse marker shapes distinguish animals within each displayed group. Display
jitter is deterministic and affects only horizontal positions. Heatmaps use
nearest-neighbor rendering with no smoothing. Summaries always use mice for
inference; there is no statistical annotation treating unit or trial counts as n.
