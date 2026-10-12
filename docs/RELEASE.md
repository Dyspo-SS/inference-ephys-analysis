# Preparing the manuscript code release

## Release scope

This repository starts with accepted unit spike counts. It includes the numerical
analysis and figure workflows, not raw-data preprocessing, sorting software,
anatomical images, or the non-electrophysiological experiments. The behavioral
Fig. 2A panel is an optional independent input. No records or numerical manuscript
results are distributed here.

## Verification status

Source syntax and repository contents can be inspected without experimental
records. This distribution has not been run end to end on the final manuscript
datasets. It must not be described as having reproduced final manuscript figures
until that run is complete. Dependencies are compatibility ranges, not a claimed
validated lockfile. The latest analysis specifications were translated into
modular code; exact numerical equivalence to earlier ad hoc script outputs is
not assumed.

Before using a release in the paper, run all three selected inputs, inspect their
completion states, exclusions, sample counts, figures and statistical tables,
and reconcile them with Methods and legends. Freeze the environment and preserve
configuration files and input hashes. Any change to inclusion, timing, model,
or correction families belongs in the release history and Methods.

## Metadata supplied by the authors

Choose the repository owner/name and an open-source license after confirming
ownership. Add actual contributor names, manuscript title/DOI, contact, and
approved data accession. A `CITATION.cff` can then be created with verified fields;
no invented author or DOI is supplied in this package.

## Upload

The source directory is ready for a repository; no remote is configured and no
upload is performed by this package. For example, from the repository root:

```bash
git init
git add .
git commit -m "Add electrophysiology analysis workflows"
# Configure the authors' actual GitHub remote, then push according to its branch policy.
```

Experimental records and outputs are ignored by default. Upload approved data
separately or explicitly revise that policy; do not accidentally commit local
recordings, patient/animal linkage files, environment secrets, or unrelated
manuscript drafts. After verification, create a versioned release and, if desired,
archive it with a persistent DOI. Cite the exact release used for the paper.

## Reporting boundaries

The code performs statistical estimation without prescribing a result direction.
Preserve nonsignificant and undefined results and report effective sample sizes.
Similarity and event locking do not establish a particular memory-retrieval
algorithm, interregional synchrony, or mediation of behavior by neural activity.
Independent group comparisons remain independent even when group-local mouse
labels happen to match.
