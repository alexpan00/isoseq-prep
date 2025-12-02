# IsoSeq Prep (p-hacking)

Are you tired of trying to use public PacBio data and not being able to run the `isoseq` pipeline because it is FASTQ or missing some tags? Here is your solution!

This repository contains a suite of tools (`p-hacking`) designed to massage, fix, and prepare PacBio sequencing data, especially data downloaded from public archives like SRA or ENA, so it plays nice with standard PacBio workflows like `isoseq refine` or `isoseq cluster`.

## Why this exists

Public data is great, but it often comes stripped of the rich metadata (BAM tags) that PacBio tools rely on. If you've ever stared at an error message complaining about missing `zm`, `rq`, or `np` tags, or wondered why your "CCS" reads are being ignored because they're in a FASTQ file, this toolkit is for you.

We provide utilities to:
- Convert FASTQ to BAM while injecting the necessary tags (`zm`, `rq`, `np`, `RG`).
- Detect and extract Kinnex concatemer signatures.
- Infer multiplexed primer sequences from data.
- Fix or add barcode tags (`bx`) to make downstream tools happy.

## Installation

The tools are packaged as a Python package. You can install it in editable mode so you can tweak things if needed:

```bash
# Clone the repo
git clone https://github.com/alexpan00/isoseq-prep.git
cd isoseq-prep

# Install dependencies and the package
pip install -e .
```

## Usage

Once installed, you have a single command `p-hacking` that exposes all the utilities.

```bash
# See all available commands
p-hacking --help

# Example: Fix a BAM file by adding ZMW and Read Quality tags
p-hacking add-zm-tag input.bam -o ready_for_isoseq.bam --compute_quality

# Example: Detect Kinnex concatemers
p-hacking detect-kinnex reads.bam --output linkers.fasta
```

## What's inside?

- **`add-zm-tag`**: The MVP. Takes a BAM (or FASTQ!) and adds the `zm` (ZMW), `rq` (Read Quality), and `np` (Num Passes) tags. It also fixes the `@RG` header so `isoseq` knows the data is CCS.
- **`detect-kinnex`**: Finds and extracts Kinnex linker sequences from concatemer reads.
- **`infer-adapters`**: Guesses primer sequences from your data if you don't have the barcode file.
- **`add-bx-tag`**: Adds barcode tags if you need to fake them or fix them.
- **`build-adapters-fasta`**: Reconstructs a primer FASTA file from BAM headers.

## Contributing

Found a new way public data is broken? Feel free to open an issue or submit a PR with a new script!

