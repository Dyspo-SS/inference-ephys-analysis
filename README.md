# Electrophysiological analyses of inferential fear

Python workflows for the electrophysiological components of Fig. 2, Fig. 4,
and V1, including animal-level statistics, plotting data, and supplementary
analyses. Inputs are event-aligned spike counts exported after spike sorting
and unit selection. The software does not perform acquisition, spike sorting,
histological localization, or LFP analysis.

## Installation

Use Python 3.10 or later in an isolated environment:

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
python -m pip install -e .
```

Dependencies are declared in `pyproject.toml`. Runtime versions, configuration,
input SHA-256, and output hashes are captured for each analysis. After a verified
manuscript run, preserve its environment with `python -m pip freeze > environment.txt`.

## Inputs

Copy the column headers from `templates/spike_counts.csv` and populate a separate
file for each dataset. Do not copy derived statistics into a count column.

| Dataset | Default input path | Configuration |
|---|---|---|
| Fig. 2 electrophysiology | `data/fig2_spike_counts.csv` | `configs/fig2.json` |
| Fig. 4 electrophysiology | `data/fig4_spike_counts.csv` | `configs/fig4.json` |
| V1 electrophysiology | `data/v1_spike_counts.csv` | `configs/v1.json` |

Each row represents **one unit, in one session, on one trial, in one time bin**.
Counts are nonnegative integers; retain valid zero-count bins. Missing recordings
are not zero counts. Units are session-specific, even if labels repeat on later
days. Animal and unit counts are read from the input, not fixed by the program.

`configs/fig2.json` selects Pairing, Conditioning, and Test for Inference, and
Conditioning and Test for No_pair. Records outside this selected stage design are
listed in `audit/exclusions.csv`. Fig. 4 selects all three stages for Control and
Gi; its GLM uses Test only. V1 uses the final Pairing session of one group and
does not run GLM or RSA. Event intervals in the configuration must match the
acquired event records. **Review the configuration before running the analysis.**

See [input schema](docs/INPUT_SCHEMA.md) for columns, identifiers, missing-data
rules, and the optional behavioral and recording-metadata files.

## Run

```bash
ephys-analysis validate --config configs/fig2.json --out results/fig2_input_check
ephys-analysis run --config configs/fig2.json --out results/fig2
ephys-analysis run --config configs/fig4.json --out results/fig4
ephys-analysis run --config configs/v1.json --out results/v1
```

Each validation or analysis run requires a new, empty output directory. A failed
run writes its error and completion status to `audit/`; it must not be presented
as a complete analysis. Use `--no-plots` to compute tables without figures.

To redraw completed results with the same configuration:

```bash
ephys-analysis plot --config configs/fig2.json --out results/fig2
```

For the optional Fig. 2A behavioral panel, populate `templates/behavior.csv` and
set `behavior.input` in `configs/fig2.json`. Without it, the exported main figure
contains panels B–D only. Behavioral animals are not matched to electrophysiology
animals. Fig. 4 output contains the electrophysiological J–L panels only.

## Analysis and outputs

- Trial-averaged unit PSTHs, animal-weighted population PSTHs, and baseline z scores.
- Animal-level window rates, exact sign-flip and group-label tests, confidence
  intervals, effect sizes, and explicit multiple-comparison families.
- Cross-validated ridge Poisson GLMs, held-out predictive gain, drop-one event
  contributions, full/reduced predictions, coefficients, fit diagnostics, and
  circular-shift null distributions.
- Population temporal-profile RSA, odd/even reliability, prespecified sensitivity
  variants, and intervention-specific unit subsampling/count-thinning controls.
- Optional independent-group behavioral contrasts, without neural–behavioral
  correlations.
- Editable SVG, vector PDF, and PNG figures; underlying CSV files are saved
  separately from statistical results.

See [analysis definitions](docs/ANALYSIS_METHODS.md),
[output inventory](docs/OUTPUTS.md), and [figure legends](docs/FIGURE_LEGENDS.md).

The analysis starts with accepted sorted units. Count-level QC cannot establish
waveform isolation or recording-location accuracy. Temporal-profile RSA compares
time-course shape; it does not track a cell ensemble across days or identify
specific retrieved memory content. Statistical conclusions must follow the
computed results, including null findings and reliability limitations.

## Repository layout

```text
configs/                Study designs and numerical analysis settings
data/                   Local input location, ignored by Git
docs/                   Methods, schemas, outputs, and release instructions
src/inference_ephys/    Analysis and plotting modules
templates/              Empty CSV headers for data preparation
pyproject.toml          Package metadata and command-line entry point
requirements.txt        Dependency ranges
```

## Publication and reuse

This distribution contains analysis code and empty schemas, not experimental
records or precomputed manuscript results. Review the
[release instructions](docs/RELEASE.md) before creating the paper's archived
release. Add the actual authors, manuscript citation, data accession, and chosen
license before public release. No license or authorship has been assigned on
behalf of the investigators.
