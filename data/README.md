# Input data

Put study CSV files in this directory or change the paths in the configuration files.
CSV paths are resolved relative to the configuration file. Recordings are not bundled.
The CSV templates in `templates/` contain column headers only.

Use deidentified animal and session identifiers. Keep all valid zero-count bins.
The pipeline never fills a missing observation with zero. `Group` is part of the
animal identifier namespace, and no pairing is inferred across independent groups.

Input files and generated results are ignored by Git by default. If an approved
data repository holds the records, add its persistent accession to the project
README before publication. A raw acquisition archive requires separate channel,
reference, event-clock, and spike-sorting metadata; binned counts alone cannot
reconstruct those records.
