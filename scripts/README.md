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
  - Notes & dependencies:
    - Requires `pysam` to read/write BAM files. If you provide FASTQ input
      the script will call `samtools import` to make a temporary BAM, so
      `samtools` must be available in your PATH when using FASTQ inputs.
    - The script will add or update an `rq` tag (read quality) if missing. By
      default it writes a constant value (configurable with
      `--default-quality`), or if `--compute_quality` is supplied the script
      computes a basewise accuracy from per-base Phred scores and writes the
      result as `rq:f:<value>`.
    - When updating/creating `@RG` entries the script ensures `DS` contains
      `READTYPE=CCS` and fills missing `PU` values so downstream Iso-Seq
      tools (e.g., `isoseq refine`) receive the expected metadata.

- `add_bx_tag.py`
  - I don't recommend anyone to use this script. I just created it because I sow
    this issue on [GitHub](https://github.com/PacificBiosciences/pbbioconda/issues/150)
    and thought that you should be able to run the tool in your data. After
    that I realize why it was a bad idea.
  - Purpose: Add or replace a per-read `bx` array tag (`bx:B:i,<len1>,<len2>`) on every
    alignment. Defaults to barcode length values 40 and 39, with CLI flags to
    override them.

  - Usage examples:

    ```bash
    # add default bx lengths
    python scripts/add_bx_tag.py input.bam

    # customize lengths and overwrite existing output
    python scripts/add_bx_tag.py input.bam --bc1-len 42 --bc2-len 38 -o with_bx.bam --force
    ```

  - Notes & dependencies:
    - Requires `pysam`.
    - Does not modify header contents; only per-read tags are updated.

- `build_adapters_fasta.py`
  - Purpose: Parse FL BAM headers to reconstruct a multiplex primer FASTA.
    Extracts barcode counts, primer indices, and sequences from read-group
    fields, fills missing indices with 25 `N`s, and writes a combined FASTA.

  - Usage example:

    ```bash
    python scripts/build_primers_fasta.py fl1.bam fl2.bam -o primers.fasta
    ```

  - Notes & dependencies:
    - Requires `pysam` to read BAM headers.
    - Validates that all inputs agree on `BarcodeCount` and do not reuse primer
      indices with conflicting sequences.

- `infer_adapters.py`
  - Purpose: Sample reads from a BAM or FASTQ file, detect polyA/polyT signals,
    and infer Iso-Seq adapter or multiplexed primer cores. With the
    `--multiplexed` flag it separates shared primer cores and tries to infer
    sample barcodes.

  - Usage examples:

    ```bash
    # infer standard adapter cores
    python scripts/infer_adapters.py sample.bam --sample 5000

    # include multiplex mode and show top 10 barcodes per end
    python scripts/infer_adapters.py sample.bam --multiplexed --top 10
    ```

  - Notes & dependencies:
    - Requires `pysam` for BAM input; FASTQ input (plain or gzipped) works
      without `pysam`.
    - Sampling and polyA/polyT heuristics are tunable via CLI flags such as
      `--sample`, `--min-poly`, and `--search-window`.
    - When multiplex mode is disabled, barcode reporting is skipped and only
      primer cores are emitted.

- `inspect_adapters.py`
  - Purpose: Highlight residual primer sequences and polyA/polyT stretches in
    individual reads pulled from a BAM file. Useful for visually verifying
    adapter trimming or spotting concatemers.

  - Usage example:

    ```bash
    python scripts/inspect_adapters.py --bam mybam.bam --max-reads 5
    ```

  - Notes & dependencies:
    - Requires `pysam` to stream reads and `edlib` for approximate matching.
    - Emits ANSI-colored output in the terminal; set a wider `--wrap-width`
      to inspect long reads.
    - Primer FASTA entries are automatically augmented with reverse
      complements before matching.

- `infer_adapters_from_concatemers.py`
  - Purpose: When you want to run `isoseq refine` but you don't have primers fasta. 
    Recover the shared primer core that persists in concatemer reads by
    scanning for terminal polyA/polyT boundaries, extracting the flanking
    adapter sequence, and building a weighted consensus.

  - Usage example:

    ```bash
    python scripts/infer_adapters_from_concatemers.py mybam.bam --sample 500000
    ```

  - Notes & dependencies:
    - Accepts BAM (requires `pysam`) or FASTQ input to supply reads.
    - Consensus heuristics are configurable via flags such as `--min-support`
      and `--n-threshold`.
    - Writes a two-entry FASTA containing the inferred 5' and 3' primer cores
      unless `--no-fasta` is provided.
    - Concatemer reads aren't really abundant so a big `--sample` is recommended

- `detect_kinnex.py`
  - Purpose: Detect Kinnex-style concatemer reads by counting multiple polyA tracts.
    It infers if the data contains Kinnex concatemers and extracts the Kinnex
    linker sequences into a FASTA file.

  - Usage example:

    ```bash
    python scripts/detect_kinnex.py input.bam --sample 1000 --output linkers.fasta
    ```

  - Notes & dependencies:
    - Requires `pysam` for BAM input and `edlib` for alignment.
    - Detects polyA hits and uses them to identify the segmentation of concatemers.
    - Outputs a FASTA file containing the inferred linker sequences.

