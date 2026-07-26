#!/usr/bin/env python
"""
extract_matrices.py

For each CRAM file and each CpG position in an input CSV, build a ~200bp
window (CpG position +/- flank bp), extract reads, build a methylation
matrix, and save one tab-delimited .bed file PER CpG position, named:
    <SAMPLE_ID>_<Chromosome>_<CpG_START>.bed
(Chromosome/CpG_START refer to the original, pre-flank CSV coordinates.)

Each CpG's output file is written to disk as soon as that CpG finishes
processing (streamed), rather than holding all matrices in memory and
writing them out at the end. This keeps memory usage low even with
thousands of CpGs in the input CSV.

Works with both BAM and CRAM files. If any CRAM files are being processed,
a reference genome FASTA must be supplied via --reference (CRAM decoding
requires it to reconstruct read sequences); this is passed through to
BamFileReadParser. BAM-only runs do not need --reference at all.
"""

import warnings

warnings.filterwarnings(
    "ignore",
    message=".*Downcasting behavior in `replace` is deprecated.*",
    category=FutureWarning,
)

import os
import sys
import glob
import logging
import argparse
import multiprocessing
from functools import partial

import pandas as pd
import numpy as np
from clubcpg.ParseBamNew import BamFileReadParser
from multiprocessing import Pool
from collections import defaultdict
import time
from pandas.core.indexes.base import InvalidIndexError

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)


def calculate_bin_coverage(input_bam, chrom, start, end, output_dir, reference_genome=None):
    """
    Take a single bin, return a matrix. This is passed to a multiprocessing Pool.
    """
    bins_no_reads = 0
    # Get reads from bam/cram file. reference_genome is required for CRAM
    # inputs (validated upstream in __main__ / process_sample) and is
    # forwarded straight through to BamFileReadParser, which uses it to
    # open the CRAM via pysam's reference_filename=.
    try:
        parser = BamFileReadParser(
            input_bam, 10, None, None, None, None, True,
            min_base_quality=5,
            reference_genome=reference_genome,
        )
    except (ValueError, FileNotFoundError) as e:
        logging.error("Could not open %s: %s" % (input_bam, e))
        return None

    try:
        reads = parser.parse_reads(chrom, start, end)
    except BaseException as e:
        # No reads are within this window, do nothing
        bins_no_reads += 1
        return None
    except:
        logging.error("Unknown error: %s:%s-%s" % (chrom, start, end))
        return None

    # Always apply correct_cpg_positions() BEFORE building the matrix, not
    # only as a fallback when the matrix comes back empty. Per the parser's
    # own docs, Bismark can report a CpG site off by 1bp even after strand
    # correction is applied -- this collapses any such (N, N+1) position
    # pairs down to N so they merge into a single matrix column instead of
    # being spread across two adjacent columns (which silently splits real
    # coverage for a CpG in half or worse).
    reads = parser.correct_cpg_positions(reads)

    try:
        matrix = parser.create_matrix(reads)
    except InvalidIndexError as e:
        logging.error("Invalid Index error when creating matrices at bin %s:%s-%s" % (chrom, start, end))
        logging.debug(str(e))
        return None
    except ValueError as e:
        logging.error("Matrix concat error at bin %s:%s-%s" % (chrom, start, end))
        logging.debug(str(e))
        return None

    # drop rows of ALL NaN
    matrix = matrix.dropna(how="all")

    if len(matrix) == 0:
        logging.info("No matrix produced at bin %s:%s-%s (no CpG-covering reads after filtering)"
                     % (chrom, start, end))

    return "%s:%s-%s" % (chrom, start, end), matrix


def _bin_worker(bin_tuple, input_bam, output_dir, sample_id, reference_genome):
    """
    Process a single CpG bin and write its output file immediately upon
    completion (runs inside a worker process). Returns a small status tuple
    instead of the matrix itself, so large matrices never need to be shipped
    back to the main process or held in memory beyond this one bin.
    """
    chrom, start, end, cpg_chrom, cpg_start = bin_tuple

    result = calculate_bin_coverage(input_bam, chrom, start, end, output_dir, reference_genome)
    if result is None:
        return (cpg_chrom, cpg_start, False)

    _, matrix = result
    if matrix is None or len(matrix) == 0:
        return (cpg_chrom, cpg_start, False)

    out_name = "%s_%s_%s.bed" % (sample_id, cpg_chrom, cpg_start)
    out_path = os.path.join(output_dir, out_name)
    matrix.to_csv(out_path, sep="\t", header=True, index=True)

    return (cpg_chrom, cpg_start, True)


def build_bins_from_csv(csv_path, flank):
    """
    Read the CpG CSV file and build (chrom, start, end) bins, centered on
    each CpG position +/- flank bp.
    """
    cpg_df = pd.read_csv(csv_path, sep=None, engine="python")  # auto-detect delimiter

    required_cols = {"chr", "start", "end"}
    missing = required_cols - set(cpg_df.columns)
    if missing:
        raise ValueError("CSV is missing required columns: %s" % missing)

    bins = []
    for _, row in cpg_df.iterrows():
        chrom = row["chr"]
        # start == end for a single CpG coordinate; use it as the center
        cpg_start = int(row["start"])
        bin_start = max(0, cpg_start - flank)
        bin_end = cpg_start + flank
        # Keep the original (pre-flank) CpG chrom/start alongside the
        # expanded bin coordinates, since output filenames are based on
        # the original CpG position, not the flanked window.
        bins.append((chrom, bin_start, bin_end, chrom, cpg_start))

    return bins


def process_sample(cram_file, bins, output_dir, threads, reference_genome=None):
    """
    Run calculate_bin_coverage over all bins for a single CRAM/BAM file.
    Each bin's output file (<SAMPLE_ID>_<Chromosome>_<CpG_START>.bed) is
    written to disk by the worker as soon as that CpG finishes -- this
    function just streams through completions and tracks progress/logging.
    """
    sample_id = os.path.splitext(os.path.basename(cram_file))[0]
    logging.info("Processing sample: %s (%d CpG bins)" % (sample_id, len(bins)))

    worker_fn = partial(
        _bin_worker,
        input_bam=cram_file,
        output_dir=output_dir,
        sample_id=sample_id,
        reference_genome=reference_genome,
    )

    n_written = 0
    n_done = 0
    if threads > 1:
        with multiprocessing.Pool(processes=threads) as pool:
            # imap_unordered streams results back as each bin finishes
            # (order doesn't matter since each bin writes its own file),
            # so files start appearing on disk immediately rather than
            # waiting for the whole sample to complete.
            for cpg_chrom, cpg_start, success in pool.imap_unordered(worker_fn, bins):
                n_done += 1
                if success:
                    n_written += 1
                else:
                    logging.debug("No output for %s:%s (sample %s)" % (cpg_chrom, cpg_start, sample_id))
                if n_done % 100 == 0:
                    logging.info("Sample %s: %d / %d bins processed, %d files written so far"
                                 % (sample_id, n_done, len(bins), n_written))
    else:
        for b in bins:
            cpg_chrom, cpg_start, success = worker_fn(b)
            n_done += 1
            if success:
                n_written += 1
            else:
                logging.debug("No output for %s:%s (sample %s)" % (cpg_chrom, cpg_start, sample_id))
            if n_done % 100 == 0:
                logging.info("Sample %s: %d / %d bins processed, %d files written so far"
                             % (sample_id, n_done, len(bins), n_written))

    logging.info("Sample %s: wrote %d / %d CpG matrix files." % (sample_id, n_written, len(bins)))


def find_cram_files(cram_dir=None, cram_file=None):
    if cram_file:
        return [cram_file]
    if cram_dir:
        files = sorted(glob.glob(os.path.join(cram_dir, "*.cram")))
        if not files:
            # also allow bam files in the same directory
            files = sorted(glob.glob(os.path.join(cram_dir, "*.bam")))
        return files
    return []


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description="Compute per-CpG-window methylation matrices from CRAM/BAM files."
    )
    arg_parser.add_argument(
        "--csv", required=True,
        help="CSV/TSV file with columns: chr, start, end, .."
    )
    arg_parser.add_argument(
        "--cram_dir", default=None,
        help="Directory containing <SAMPLE_ID>.cram (or .bam) files to process"
    )
    arg_parser.add_argument(
        "--cram", default=None,
        help="Path to a single CRAM/BAM file (alternative to --cram_dir)"
    )
    arg_parser.add_argument(
        "--reference", default=None,
        help="Reference genome FASTA (must have a .fai index alongside it). "
             "REQUIRED if any of the input files are CRAM -- CRAM decoding "
             "needs the reference to reconstruct read sequences. Not needed "
             "for BAM-only input."
    )
    arg_parser.add_argument(
        "-o", "--output", default=None,
        help="Folder to save output matrices (default: same dir as input CRAM)"
    )
    arg_parser.add_argument(
        "--flank", type=int, default=100,
        help="Number of bp to add upstream/downstream of each CpG position (default: 100)"
    )
    arg_parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of worker processes for parallel bin processing (default: 4)"
    )
    args = arg_parser.parse_args()

    if not args.cram_dir and not args.cram:
        arg_parser.error("You must provide either --cram_dir or --cram")

    cram_files = find_cram_files(cram_dir=args.cram_dir, cram_file=args.cram)
    if not cram_files:
        logging.error("No CRAM/BAM files found.")
        sys.exit(1)

    # If any of the input files are CRAM, a reference genome is required.
    # Fail fast with a clear message rather than letting BamFileReadParser
    # raise deep inside a worker process for every single bin.
    cram_present = any(f.endswith(".cram") for f in cram_files)
    if cram_present and not args.reference:
        arg_parser.error(
            "One or more input files are CRAM (%s), which requires a "
            "reference genome to decode. Please supply --reference "
            "/path/to/genome.fa (with a .fai index alongside it)."
            % ", ".join(f for f in cram_files if f.endswith(".cram"))
        )
    if args.reference and not os.path.exists(args.reference):
        arg_parser.error("Reference genome file not found: %s" % args.reference)

    # Set output dir
    if not args.output:
        output_folder = os.path.dirname(cram_files[0]) or "."
    else:
        output_folder = args.output

    try:
        os.makedirs(output_folder, exist_ok=True)
    except FileExistsError:
        pass

    logging.info("Building bins from CSV: %s (flank=%dbp)" % (args.csv, args.flank))
    bins = build_bins_from_csv(args.csv, args.flank)
    logging.info("Total bins to process per sample: %d" % len(bins))

    for cram_file in cram_files:
        process_sample(
            cram_file=cram_file,
            bins=bins,
            output_dir=output_folder,
            threads=args.threads,
            reference_genome=args.reference,
        )

    logging.info("Done.")