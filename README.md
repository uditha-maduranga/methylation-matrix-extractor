# methylation-matrix-extractor

Fast and flexible extraction of read-level DNA methylation matrices from BAM/CRAM files for arbitrary genomic regions.

Given a list of CpG coordinates and one or more single-ended BAM/CRAM files, this tool builds a small window around each CpG, pulls the overlapping reads, and writes out a read-by-CpG methylation matrix as a tab-delimited `.bed` file — one file per CpG position, per sample. Output is streamed to disk as each CpG finishes, so memory usage stays flat even with CSV inputs containing thousands of CpG sites.

---

## Table of contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input format](#input-format)
- [Usage](#usage)
- [Command-line arguments](#command-line-arguments)
- [Output format](#output-format)
- [How it works](#how-it-works)
- [Performance notes](#performance-notes)
- [Known limitations](#known-limitations)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Features

- Works with **both BAM and CRAM** input files.
- Accepts **single or multiple samples** (single file, or a whole directory).
- Builds a configurable window (`± flank` bp) around each CpG position.
- Streams output to disk **per CpG, as soon as it's ready** — no need to hold thousands of matrices in memory.
- **Parallelized** across CpG bins via `multiprocessing`.
- Merges methylation information of reads from forward and reverse strands and shows together in the matrix.
- Simple, predictable output naming: `<SAMPLE_ID>_<Chromosome>_<CpG_START>.bed`.

## Requirements

- Linux (tested environment; other Unix-likes likely work)
- Conda / Miniconda / Mambaforge
- Python 3.9 (pinned in the provided environment)

Core dependencies (see [`clubcpg_environment.yml`](./clubcpg_environment.yml) for exact pinned versions):

| Package | Purpose |
|---|---|
| [`clubcpg`](https://pypi.org/project/clubcpg/) | Provides the base package structure — see the important note below, though |
| `pysam` | BAM/CRAM I/O |
| `pandas` / `numpy` | CSV parsing and matrix handling |
| `scikit-learn`, `scipy`, `fastcluster`, `joblib`, `pebble`, `tqdm` | Pulled in as `clubcpg` dependencies |

> **Important:** this repo does **not** use the stock parser that ships with the `clubcpg` PyPI package (`clubcpg/ParseBam.py`). It uses a custom, modified parser — `ParseBamNew.py`, included in this repository under [`clubcpg_modification/ParseBamNew.py`](./clubcpg_modification/ParseBamNew.py) — which `extract_matrices.py` imports as `clubcpg.ParseBamNew`. This custom parser adds CRAM support, single-ended read compatibility, and forward/reverse-strand methylation merging on top of the stock module. **It must be manually installed into the `clubcpg` package directory after creating the conda environment** — see [Installation](#installation) step 4 below. This is not optional; the script will fail to import without it.

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/uditha-maduranga/methylation-matrix-extractor.git
   cd methylation-matrix-extractor
   ```

2. Create the conda environment from the provided YAML file:

   ```bash
   conda env create -f clubcpg_environment.yml
   ```

3. Activate it:

   ```bash
   conda activate clubcpg_env
   ```

4. **Add the custom modified parser (`ParseBamNew.py`).** The `clubcpg` package installed via pip only ships its default `ParseBam.py`; the modified `ParseBamNew.py` in this repo (`./clubcpg_modification/ParseBamNew.py`) needs to be added into that same installed package directory, since that's where `extract_matrices.py` imports it from (`from clubcpg.ParseBamNew import BamFileReadParser`).

   First, find where `clubcpg` is installed inside your active environment:

   ```bash
   python -c "import clubcpg, os; print(os.path.dirname(clubcpg.__file__))"
   ```

   This will print something like:

   ```
   /cluster/home/<YOUR_FOLDER>/anaconda3/envs/clubcpg_env/lib/python3.9/site-packages/clubcpg
   ```

   Then move `ParseBamNew.py` into that directory using `mv`:

   ```bash
   mv ./clubcpg_modification/ParseBamNew.py \
      /cluster/home/<YOUR_FOLDER>/anaconda3/envs/clubcpg_env/lib/python3.9/site-packages/clubcpg/ParseBamNew.py
   ```

   (Replace the destination path with whatever the `python -c ...` command above actually printed for your system — the `<YOUR_FOLDER>` and Python version segment will differ depending on your username and environment location.)

   > This step needs to be repeated any time you recreate the `clubcpg_env` environment from scratch, since a fresh `conda env create` only reinstalls the stock `clubcpg` package without `ParseBamNew.py`.

## Input format

A CSV (or TSV — the delimiter is auto-detected) listing one CpG position per row. **`start` and `end` must be equal**, representing the single-base cytosine in the CpG coordinate:

```
chr,start,end
chr1,136522,136522
chr1,978101,978101
```

Required columns: **`chr`, `start`, `end`**. Any additional columns are ignored by the script and can be left in the file if convenient (e.g. if the CSV originates from a differential methylation analysis with extra annotation columns).

## Usage

Process every CRAM/BAM file in a directory against a CpG list. **A reference genome FASTA (with a `.fai` index alongside it) is required whenever any CRAM files are involved:**

```bash
python extract_matrices.py \
  --csv cpg_sites.csv \
  --cram_dir /path/to/cram_files/ \
  --reference /path/to/genome.fa \
  --output /path/to/output/ \
  --flank 100 \
  --threads 8
```

Process a single sample:

```bash
python extract_matrices.py \
  --csv cpg_sites.csv \
  --cram /path/to/sample1.cram \
  --reference /path/to/genome.fa \
  --output /path/to/output/ \
  --flank 100 \
  --threads 8
```

If you're processing **BAM files only**, `--reference` can be omitted entirely:

```bash
python extract_matrices.py \
  --csv cpg_sites.csv \
  --cram_dir /path/to/bam_files/ \
  --output /path/to/output/ \
  --flank 100 \
  --threads 8
```

## Command-line arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--csv` | Yes | — | CSV/TSV file with (at minimum) `chr`, `start`, `end` columns |
| `--cram_dir` | One of `--cram_dir` / `--cram` required | — | Directory containing `<SAMPLE_ID>.cram` (or `.bam`) files. All `*.cram` files are used if present; `*.bam` files are only used as a fallback if no CRAMs are found in the directory |
| `--cram` | One of `--cram_dir` / `--cram` required | — | Path to a single CRAM/BAM file, as an alternative to `--cram_dir` |
| `--reference` | **Required if any input file is CRAM** | `None` | Reference genome FASTA, with a `.fai` index alongside it (e.g. generated via `samtools faidx genome.fa`). Needed to decode CRAM reads back into sequences. Not required for BAM-only input — the script checks upfront and exits with a clear error if a CRAM is present without `--reference` |
| `-o`, `--output` | No | Same directory as the input CRAM/BAM | Folder to write output `.bed` files to (created if it doesn't exist) |
| `--flank` | No | `100` | Number of bp added upstream and downstream of each CpG position, giving a total window of `2 × flank` bp |
| `--threads` | No | `4` | Number of worker processes used to parallelize across CpG bins. Set to `1` to run single-threaded/sequentially |

## Output format

For every sample and every CpG position that produces at least one non-empty matrix row, one file is written:

```
<SAMPLE_ID>_<Chromosome>_<CpG_START>.bed
```

- `SAMPLE_ID` is the CRAM/BAM filename without its extension.
- `Chromosome` and `CpG_START` refer to the **original CSV coordinates** (not the flanked window).

Each file is a tab-delimited table:
- **Rows**: individual reads covering the window.
- **Columns**: CpG positions within the window.
- **Values**: methylation call at each CpG for each read (matrix format produced by `BamFileReadParser.create_matrix`).

Rows that are entirely missing (no CpG in the window is covered by that read) are dropped before writing. If a CpG window has no usable reads at all, **no file is written** for that CpG/sample — this is expected behavior for low-coverage regions, not an error.

## How it works

1. **Build bins** (`build_bins_from_csv`): each CpG's `start` position becomes the center of a window `[start − flank, start + flank]` (clamped at 0).
2. **Per-bin extraction** (`calculate_bin_coverage`, run in a worker process):
   - Opens the CRAM/BAM via `BamFileReadParser`, passing along `--reference` when supplied (required and validated upfront for CRAM input; ignored for BAM).
   - Parses reads overlapping the window.
   - Applies `correct_cpg_positions()` **before** matrix construction, collapsing any adjacent (N, N+1) CpG position artifacts into a single column so coverage for one real CpG isn't split across two.
   - Builds the read × CpG methylation matrix and drops all-NaN rows.
3. **Streamed write** (`_bin_worker`): as soon as a bin's matrix is ready, it's written straight to disk as its own `.bed` file — matrices are never accumulated across bins.
4. **Parallel dispatch** (`process_sample`): bins are distributed across worker processes with `multiprocessing.Pool.imap_unordered`, so files appear on disk as soon as each CpG finishes, in whatever order completes first. Progress is logged every 100 processed bins.
5. This repeats independently for every sample found via `--cram_dir` / `--cram`.

## Performance notes

- Each CpG bin currently opens its own `BamFileReadParser` instance (i.e. the CRAM/BAM file handle is opened per bin, not once per worker). For very large CpG lists this reopening overhead can dominate runtime — if extraction feels slow, this is the first place to optimize (e.g. via a `Pool` initializer that opens the parser once per worker process).
- `--threads` controls parallelism across **bins within one sample**; samples themselves are still processed one after another.
- Runtime scales with the number of CpGs × samples × sequencing depth in each window: wider `--flank` values or higher-coverage samples will increase per-bin processing time.

## Known limitations

- **No resume support.** Re-running the script will reprocess every CpG and overwrite any existing `.bed` files; there's no automatic skip-if-exists check.
- **CRAM/BAM mixing:** if a `--cram_dir` contains both CRAM and BAM files, only the CRAM files are processed (BAM is used only as a directory-wide fallback when zero CRAMs are found).
- Designed and tested for **single-ended** data.
- Errors during read parsing (e.g. genuinely malformed regions) and simply-empty windows (e.g. low coverage) are currently logged similarly, so it's worth checking `DEBUG`-level logs if you suspect a specific region is failing rather than just lacking coverage.

## Troubleshooting

**`ModuleNotFoundError: No module named 'clubcpg.ParseBamNew'` (or similar import error):**
This means step 4 of [Installation](#installation) hasn't been done (or was undone by recreating the environment). Locate your `clubcpg` install directory with `python -c "import clubcpg, os; print(os.path.dirname(clubcpg.__file__))"` and copy `clubcpg_modification/ParseBamNew.py` into it.

**Script exits immediately with an error about `--reference` being required:**
This is expected and intentional — the script checks upfront whether any input file is CRAM and requires `--reference` in that case, rather than failing deep inside a worker process on every single bin. Supply `--reference /path/to/genome.fa`, making sure a matching `.fai` index exists alongside it (`samtools faidx genome.fa`).

**A CpG produces no output file:**
This is expected when there are no reads covering that window after filtering — check the `DEBUG`-level logs (`logging.debug`) for that chrom/position to confirm.

**`InvalidIndexError` or matrix concat errors in the log:**
These are caught and logged per-bin (with the offending `chrom:start-end`); that bin is simply skipped and processing continues with the rest.

## License

MIT