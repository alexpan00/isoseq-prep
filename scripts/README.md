# Scripts

This directory contains small helper scripts used to prepare PacBio sequencing
data (e.g., files downloaded from public repositories or sequencing providers).

Currently included scripts

- `add_zm_tag.py`
  - Purpose: Read a BAM file and add/replace a per-read integer tag `zm` (as
    `zm:i:N`) derived from the read name. The script extracts the number after
    the first `/` in the read name and writes it into the `zm` tag. It also
    ensures that the BAM header's `@RG` records include a `PU` (platform unit)
    tag — derived from the first part of the first read name (the substring
    before the first `/`). If no `@RG` lines exist, a minimal `@RG` with
    `ID:RG0` and the computed `PU` is added.

  - Why this script exists: I struggled with sample `ERR13885894` where the
    `zm` and `PU` information was missing; I created this script to
    populate the `zm` tag (and add `PU` to `@RG`) so downstream tools can use
    those fields.

  - Usage examples:

    ```bash
    # basic: create an output BAM with zm tags added
    python scripts/add_zm_tag.py /path/to/input.bam

    # specify output, show progress, overwrite existing output
    python scripts/add_zm_tag.py /path/to/input.bam -o /path/to/out.zm.bam --verbose --force
    ```

  - Notes & dependencies:
    - Requires `pysam` to read/write BAM files. `samtools` is useful for
      validating/indexing the output BAM.
    - The script sets `zm` to `0` when the read name doesn't contain a
      numeric field after the first `/`.
    - PU is derived from the first alignment's query name. If that's not the
      correct heuristic for your data, consider adjusting the script or
      providing a constant PU and/or mapping logic.

