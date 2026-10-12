# Input schema

## Required spike-count table

| Column | Type | Meaning |
|---|---|---|
| Mouse_ID | string | Animal identifier within Group; preserve leading zeros |
| Session_ID | string | Recording session identifier within animal |
| Group | string | Exact key in the configuration's groups object |
| Region | string | Exact configured region label |
| Unit_ID | string | Sorted unit identifier within region and session |
| Stage | string | Exact configured stage label |
| Trial_ID | string | Trial identifier within session |
| Trial_Order | positive integer | Original acquisition order, not a renumbered order after exclusions |
| Time_Bin_Start_s | number | Bin start relative to the first cue of the trial |
| Time_Bin_End_s | number | Bin end on the same clock |
| Exposure_s | positive number | Observed bin duration in seconds |
| Spike_Count | nonnegative integer | Number of accepted spikes in this bin, including zero |
| Sound | 0 or 1 | Bin entirely inside the actual sound interval |
| Light | 0 or 1 | Bin entirely inside the actual light interval |
| Shock | 0 or 1 | Bin entirely inside the actual shock interval |
| Valid_Bin | 0 or 1 | Whether the observation is valid after acquisition/artifact review |

The unique key is Group + Mouse_ID + Session_ID + Region + Unit_ID + Stage +
Trial_ID + Time_Bin_Start_s. No separate unit-count file is needed: unit numbers
are derived from these keys. Additional columns such as Day may be kept but do
not determine the analysis. No data-type column is required or overwritten.

Use a composite identifier if the same recording identifier was reused within
one animal. Reuse of Unit_ID across days does not establish cell identity.

## Epochs and event alignment

The provided configurations select 250-ms complete bins and a −5 to 0 s baseline.
Default intervals are:

| Stage | Epoch | Sound | Light | Shock |
|---|---|---|---|---|
| Pairing | −5 to 25 s | 0 to 10 s | 10 to 20 s | absent |
| Conditioning | −5 to 17 s | absent | 0 to 10 s | 10 to 12 s |
| Test | −5 to 15 s | 0 to 10 s | absent | absent |

These intervals are an explicit study configuration, not a way to fabricate
event labels. Recover timestamps from the acquisition event stream and verify
the configuration against them. Trials with variable cue timing require
event-specific alignment and a revised model implementation; this package
requires equal bin grids and fixed event timing within each selected stage.
Partially observed bins and partial event overlap are not supported.

Intertrial intervals are not concatenated into trial-relative time. Event gaps
and baseline periods must come from their true acquired timestamps. No counts
are generated for unrecorded intervals. Absolute event times should remain in
the upstream recording archive even though they are not columns in this table.

## Exclusions and missing records

Default behavior is to stop on invalid/incomplete unit-trial epochs. The
alternative `drop_unit_trial` policy removes the entire affected unit-trial,
records the reason, and preserves the original Trial_Order. It must be chosen
before evaluating results. Negative or noninteger counts and event-label
conflicts are data errors and always stop the run.

Units with zero baseline SD remain in unstandardized rate analyses but do not
contribute to z-score profiles. All-zero units or folds cannot provide the
specified GLM gain; model failures are recorded. The default GLM policy stops;
an explicit `exclude` policy permits analysis of estimable units with a full
failure table. These are separate analysis-specific inclusion rules, not a
replacement for sorting QC. Planned but unestimable tests remain in correction
families rather than reducing the multiplicity denominator.

The implementation permits different trial availability across units after
explicit exclusions, but preserves one session per mouse, region, and stage.
Multiple sessions require a prespecified session selection/aggregation method.
Two groups are treated as independent, never as paired because numeric IDs match.

## Optional behavioral table

`Mouse_ID, Group, Day, Freezing_Percent` with one session-level percentage per
animal/day. Whole-session freezing and cue-only freezing must not be mixed.
The shipped optional design compares Inference with each control on each of
eight days, using Welch tests and a single 16-test Holm family. Edit the day/group
configuration if the actual behavioral design differs. No animal linkage to
electrophysiology is inferred. The table is used only when its path is configured.

## Acquisition and sorting metadata

`recording_metadata.csv` and `unit_quality_metadata.csv` are header-only archive
templates; they are **not automatically analyzed**. Fill their fields from the
acquisition system and actual sorting QC. Do not infer channel count, ISI
violations, waveform stability, cell type, or isolation from spike-count bins.
No missing sorting-quality values should be reported as zero.

Animal, session, unit-session, trial, and channel counts are different quantities.
Session totals in the audit count shared dCA1/ACC recording sessions once.
Regional sample tables should not be summed to count sessions or animals.
