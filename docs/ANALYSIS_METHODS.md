# Analysis definitions

## Sampling and rate summaries

The independent unit is the mouse. One selected session is used per mouse,
region, and stage. No cell identity is matched across days. A bin rate is
Spike_Count / Exposure_s. Complete valid bins are retained, including zero counts.
PSTHs average trials within a unit-session, units within a mouse, then mice within
a group. Animal weighting is equal, independent of unit yield.

Baseline mean and sample SD are calculated from all retained single-trial
baseline bins of each unit-session. Baseline z = (trial-averaged bin rate − baseline
mean) / baseline SD. Zero SD yields undefined z, not an added pseudocount. Neither
temporal nor interunit smoothing is applied. Event changes use each trial's
pre-first-cue baseline, then average trials, units, and mice in that order.

Intervention primary measures are test-sound rate and the natural logarithm of
the animal sound/baseline rate ratio. The ratio is computed after unit averaging.
If either rate is zero the log ratio is undefined; no epsilon is added. Pre-cue
rate and absolute change are secondary. Both epochs are post-treatment. V1's
primary measure is sound-window minus baseline firing; early and late windows
are secondary and light responses descriptive. No duration endpoint is tested.

## GLM

Separate unit-session Poisson models use log expected count = log Exposure_s +
intercept + standardized trial order + coefficients for configured stimulus
intervals. The fixed ridge penalty is 0.5 × lambda × sum of squared nonintercept
coefficients. The reduced model omits all event terms. Drop-one contributions
omit one event at a time and are descriptive.

Complete odd trials train and complete even trials validate, then the roles
reverse. Conditioning uses leave-one-complete-trial-out and is exploratory.
Trial-order mean and SD are estimated on the training partition only. Counts
are not z scored for fitting. Folds are pooled before normalization:

`gain (bits/spike) = sum(LL_full − LL_reduced) / [sum(heldout spikes) × ln(2)]`.

Zero-spike validation folds still contribute their log-likelihood difference;
their separate bits/spike value is undefined. Unit gains are averaged within
mouse, followed by animal-level group summaries. Likelihood constants cancel
between models and are omitted in exported model scores. Fitting aggregates
bins with identical design rows; this is an exact sufficient-statistic
aggregation, not downsampling of the spike process.

For each of 999 default null draws, counts are circularly shifted within each
trial relative to the event indicators. The same draw and trial shift are shared
across units and regions recorded together. Zero shifts are allowed. Null full
models are refitted with the same folds and penalty. The reduced model can be
reused because its predictors are constant within a trial and its sufficient
statistics are unchanged by a circular shift. The upper-tail probability is
`(1 + number(null gain >= observed gain)) / (number of draws + 1)`.

This is a test of global event-time locking, conditional on the circular-shift
null, not a test that average gain is greater than zero and not a between-group
effect. Trial structure, nonstationarity, and the small conditioning trial count
limit interpretation. Predictive-gain confidence intervals and shift P values
answer different questions. No synchrony, connectivity, or memory-content
selectivity is inferred from the event model.

Optimizer convergence, training Pearson dispersion, coefficients, trial splits,
held-out predictions, and model failures are exported. A ridge-Poisson fit does
not remove biological overdispersion; group inference uses animals rather than
Poisson coefficient standard errors. The pipeline does not optimize a model
family or analysis window for significance.

## Temporal-profile RSA

For each mouse and region, independently average normalized units within each
day. The default temporal feature vector contains eight 250-ms bins over 0–2 s
after each cue onset. Pearson correlation between two vectors quantifies their
time-course similarity. PS, CL, and TS denote pairing sound, conditioning light,
and test sound. All available matrix edges are exported.

CL–TS is the primary edge. The association comparison has shared CL and TS
recordings in both groups; PS edges are descriptive and available only in the
group with selected Pairing recordings. Intervention groups have three-condition
matrices, with PS edges descriptive. Group matrices average animal r values.
Group inference uses differences of mean Fisher z = atanh(r), with epsilon
clipping recorded and reported for correlations too close to ±1. No theoretical
RDM model fitting or tracked-unit representational analysis is implemented.

Odd and even trials receive separate baseline estimates and normalization;
within-condition profile correlations describe reliability. Low reliability
is not an exclusion criterion. Differences in reliability, onset latency,
adaptation, and generic sensory/behavioral states can influence similarity.

Sensitivity variants retain 0–1, 0–2, and 0–3 s z profiles at 250 ms; a 0–2 s
profile formed by averaging adjacent standardized bins to 500 ms; and a raw-rate
0–2 s profile. The wider-bin variant does not re-estimate baseline SD.

## Intervention count controls

Sample the configured number of units without replacement, independently in CL
and TS, for each animal and region. The default is six units and 200 repetitions.
Animals lacking this number in either stage are listed as ineligible, not silently
replaced with duplicates. Additional Control test counts, including baseline,
are binomially thinned with probability min(1, Gi/Control group-mean sound-rate
ratio), computed from the actual input. Normalization and correlations are
recomputed after thinning. Group summaries use animal median correlations over
repetitions; repetition quantiles are Monte Carlo intervals, not animal CIs.
No confirmatory P values are assigned to these controls. They do not match
response reliability or establish that an effect is independent of reliability.

## Animal-level inference

Paired changes use two-sided sign flips, requiring sign exchangeability under
the null. Independent groups use two-sided group-label permutations of mean
differences. Exact enumeration is used up to the configured assignment limit;
larger designs use the configured Monte Carlo count and a plus-one correction.
The output labels the method and assignment count; it never assumes a fixed
number of animals. Independent-group exchangeability is an assumption, not
guaranteed by unequal observational group variances.

Percentile bootstrap CIs resample animals, separately within each group for
group contrasts, with 10,000 default draws. Intervals are not multiplicity
adjusted. Paired d_z is mean change / sample SD of changes. Hedges g uses pooled
SD and correction 1 − 3 / [4(n_A+n_B) − 9]. RSA's primary effect is the mean
Fisher-z difference; plotted group intervals remain on the r scale.

### Default correction families

| Dataset | Family | Default planned tests | Correction |
|---|---|---|---|
| Fig. 2 | Group × region × available stimulus response | 16 | Holm |
| Fig. 2 | Group × region × available stage GLM event locking | 10 | BH FDR |
| Fig. 2 | Regional CL–TS contrasts | 2 | Holm |
| Fig. 2 | Region × five RSA definitions | 10 | Holm |
| Fig. 4 | Sound rate and log ratio × region | 4 | Holm |
| Fig. 4 | Baseline rate and absolute change × region | 4 | Holm |
| Fig. 4 | Group × region Test GLM event locking | 4 | Holm |
| Fig. 4 | Regional CL–TS contrasts | 2 | Holm |
| Fig. 4 | RSA variants and count controls | descriptive | none |
| V1 | Full sound versus baseline | 1 | none |
| V1 | Early and late sound versus baseline | 2 | Holm |

Planned hypotheses are determined by the study configuration, never by the
observed P values. Unavailable tests retain a family slot and receive NA in
their reported result; a P value of 1 is used internally for correction only.
Changing stages, windows, or regions changes the planned family and must be
documented. No Group × Region interaction, causal mediation, or neuronal
correlation with separately identified behavioral animals is tested.
