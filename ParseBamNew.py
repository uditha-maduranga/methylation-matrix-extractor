import os
import pysam
import pandas as pd
from collections import defaultdict
import logging
import re


class BamFileReadParser:
    """
    Used to simplify the opening and reading from BAM/CRAM files. Files must be coordinate sorted and indexed.

    :Example:
        >>> from clubcpg.ParseBamNew import BamFileReadParser
        >>> parser = BamFileReadParser("/path/to/data.BAM", quality_score=20, read1_5=3, read1_3=4, read2_5=7, read2_3=1)
        >>> reads = parser.parse_reads("chr7", 10000, 101000)
        >>> reads = parser.correct_cpg_positions(reads) # This step is optional, but highly recommended
        >>> matrix = parser.create_matrix(reads)


    """

    def __init__(self, bamfile, quality_score, read1_5=None, read1_3=None,
                 read2_5=None, read2_3=None, no_overlap=True, apply_strand_shift=False,
                 min_base_quality=5, reference_genome=None):
        """
        Class used to read WGBSeq reads from a BAM/CRAM file, extract methylation, and convert into data frame

        :param bamfile: Path to bam/cram file location
        :param quality_score: Only include reads >= this MAPQ (matches MethylDackel's -q, default 10)
        :param read1_5: mbias ignore read1 5'
        :param read1_3: mbias ignore read1 3'
        :param read2_5: mbias ignore read2 5'
        :param read2_3: mbias ignore read2 3'
        :param no_overlap: bool. If overlap exists between two reads, ignore that region from read 2.
        :param apply_strand_shift: bool. If True, reads whose XG tag is 'GA'
            (i.e. matched the G->A converted / bottom-strand reference) have
            their CpG position shifted by -1 to align with the top-strand C
            coordinate. This is the classic Bismark raw-BAM convention.
            DEFAULT IS FALSE: empirical testing (diagnose_shift.py) on this
            pipeline's CRAM files showed positions are ALREADY normalized to
            the top-strand C for every read regardless of originating
            strand -- applying the shift on top of that artificially SPLITS
            correctly-merged CpG coverage into two adjacent columns instead
            of fixing anything. If you use this parser on a different
            dataset (e.g. a raw/unprocessed Bismark BAM), re-verify with
            diagnose_shift.py before assuming either setting is correct.
        :param min_base_quality: Minimum per-base Phred quality to count an
            individual CpG call (matches MethylDackel's -p, default 5).
            A base with quality below this is dropped even if the read as a
            whole passed the MAPQ filter.
        :param reference_genome: Path to the reference genome FASTA (must
            have a .fai index alongside it, e.g. via `samtools faidx`).
            REQUIRED if bamfile is a CRAM file, since CRAM decoding needs the
            reference to reconstruct read sequences. Ignored for BAM files.
        """

        self.mapping_quality = quality_score
        self.min_base_quality = min_base_quality
        self.bamfile = bamfile
        self.reference_genome = reference_genome
        self.read1_5 = read1_5
        self.read1_3 = read1_3
        self.read2_5 = read2_5
        self.read2_3 = read2_3
        self.full_reads = []
        self.read_cpgs = []
        self.no_overlap = no_overlap
        self.apply_strand_shift = apply_strand_shift

        if read1_5 or read2_5 or read1_3 or read2_3:
            self.mbias_filtering = True
        else:
            self.mbias_filtering = False

        is_cram = bamfile.endswith(".cram")

        if is_cram:
            if not reference_genome:
                raise ValueError(
                    "A reference genome FASTA is required to open CRAM file "
                    "'%s', but none was provided. Pass reference_genome="
                    "/path/to/genome.fa to BamFileReadParser (and make sure "
                    "a .fai index exists alongside it, e.g. via "
                    "`samtools faidx genome.fa`)." % bamfile
                )
            if not os.path.exists(reference_genome):
                raise FileNotFoundError(
                    "Reference genome file not found: %s" % reference_genome
                )

        # NOTE: mode 'rb' works for BAM. For CRAM, we explicitly use mode
        # 'rc' and pass reference_filename= (validated above), since CRAM
        # decoding requires the reference genome to reconstruct read
        # sequences, and relying on htslib's online reference lookup can
        # silently fail / time out in environments without internet access
        # (e.g. HPC compute nodes).
        open_kwargs = {}
        if is_cram:
            open_kwargs["reference_filename"] = reference_genome

        self.OpenBamFile = pysam.AlignmentFile(
            bamfile,
            "rc" if is_cram else "rb",
            **open_kwargs
        )
        # Check for presence of index file
        index_present = self.OpenBamFile.check_index()
        if not index_present:
            raise FileNotFoundError("BAM/CRAM file index is not found. Please create it using samtools index")


    # From open bam file, get locaiton of first read from the provided chromosome
    def get_location_of_first_read(self, chromosome):

        # Get reference lenghts
        ref_lens = dict(zip(self.OpenBamFile.references, self.OpenBamFile.lengths))

        for read in self.OpenBamFile.fetch(chromosome, 0, ref_lens[chromosome]):
            reads_start_loc = read.reference_start
            break

        return reads_start_loc

    # Get reads from the bam file, extract methylation state
    def parse_reads(self, chromosome: str, start:int , stop: int):
        """
        :param chromosome: chromosome as "chr6"
        :param start: start coordinate
        :param stop: end coordinate
        :return: List of reads and their positional tags as assigned by bismark
        """
        reads = []
        for read in self.OpenBamFile.fetch(chromosome, start, stop):
            # Filter on mapping quality AND exclude duplicate / secondary /
            # supplementary / QC-fail alignments. This matches MethylDackel's
            # default ignoreFlags=0xF00, which masks out secondary (0x100),
            # QC-fail (0x200), duplicate (0x400), and supplementary (0x800)
            # reads.
            if (read.mapping_quality >= self.mapping_quality
                    and not read.is_duplicate
                    and not read.is_secondary
                    and not read.is_supplementary
                    and not read.is_qcfail):
                reads.append(read)

        ## CIGAR FILTERING BY C. COARFA
        read_cpgs = []
        self.skipped_reads = set()

        read_index = -1
        self.query_count_hash = {}
        for read in reads:
            read_index +=1

            if not (read.query_name in self.query_count_hash):
                self.query_count_hash[read.query_name]=0
            
            self.query_count_hash[read.query_name] += 1
            # if (self.query_count_hash[read.query_name]>2):
            #     logging.info("Found read with more than 2 mappings: %s --> %s\n"%(read.query_name, self.query_count_hash[read.query_name]))

            # NOTE: We used to reject the entire read here if its CIGAR
            # contained a deletion ('D') or reference skip ('N'), because a
            # naive zip(get_aligned_pairs(), XM) breaks across those
            # operations (they consume a reference base with no
            # corresponding query base / XM character, which used to
            # misalign everything downstream). MethylDackel does NOT reject
            # such reads -- it just correctly walks past the deleted/skipped
            # bases. We now do the same: index directly into the XM string
            # and query qualities by query position (qpos) rather than
            # relying on parallel iteration order, so D/N operations are
            # naturally skipped (they yield qpos=None) without breaking
            # alignment for the rest of the read. This also lets us apply
            # MethylDackel's per-base Phred quality filter (-p, default 5),
            # which requires indexing qualities by qpos too.
            if read.cigarstring is not None:

                reduced_read = []

                # Determine whether this read's methylation calls need the
                # -1 shift to align to the top-strand C position of each CpG.
                # Gated behind self.apply_strand_shift (see __init__ docstring
                # for why this defaults to False for this pipeline's data).
                if self.apply_strand_shift:
                    try:
                        needs_shift = (read.get_tag('XG') == 'GA')
                    except KeyError:
                        # Fallback if XG tag is unavailable for some reason:
                        # use flag-based inference (directional-library
                        # assumption only).
                        if read.is_paired:
                            needs_shift = (read.is_read1 and read.is_reverse) or \
                                          (read.is_read2 and not read.is_reverse)
                        else:
                            needs_shift = read.is_reverse
                else:
                    needs_shift = False

                xm = read.get_tag('XM')
                qualities = read.query_qualities  # indexed by qpos; None if unavailable

                for qpos, rpos in read.get_aligned_pairs():
                    if qpos is None:
                        # Deletion ('D') or reference skip ('N') -- no query
                        # base / XM character exists for this position.
                        # Nothing to record; crucially, this does NOT
                        # consume an XM character, so subsequent bases stay
                        # correctly aligned.
                        continue
                    if rpos is None:
                        # Insertion ('I') or soft-clip ('S') -- no reference
                        # position to record this base against.
                        continue

                    # Per-base Phred quality filter (MethylDackel's -p,
                    # default minPhred=5). A base can fail this even if the
                    # read as a whole passed the MAPQ filter.
                    if qualities is not None and qualities[qpos] < self.min_base_quality:
                        continue

                    tag = xm[qpos]
                    if needs_shift:
                        reduced_read.append((rpos - 1, tag))
                    else:
                        reduced_read.append((rpos, tag))
    
                # if MBIAS was set, slice the joined list
                if self.mbias_filtering:
                    if read.is_read1:
                        mbias_5_prime = self.read1_5
                        # note taking the NEGATIVE of the value for the 3-prime
                        mbias_3_prime = -self.read1_3
                        if mbias_3_prime == 0:
                            mbias_3_prime = None
                        reduced_read = reduced_read[mbias_5_prime:mbias_3_prime]
                    if read.is_read2:
                        mbias_5_prime = self.read2_5
                        mbias_3_prime = -self.read2_3
                        if mbias_3_prime == 0:
                            mbias_3_prime = None
                        reduced_read = reduced_read[mbias_5_prime:mbias_3_prime]
                    
                read_cpgs.append(reduced_read)
            else:
                self.skipped_reads.add(read.query_name)

        self.full_reads = reads
        self.read_cpgs = read_cpgs

        # Correct overlapping paired reads if set, this is default behavior
        if self.no_overlap:
            try:
                read_cpgs = self.fix_read_overlap(reads, read_cpgs)
            except AttributeError:
                pass
                # print("Could not determine read 1 or 2. {}:{}-{}".format(chromosome, start, stop))
                # sys.stdout.flush()

        # Filter the list for positions between start-stop and CpG (Z/z) tags
        output = []
        read_cpg_index = -1
        
        found_cpg_count = 0
        
        for read_cpg in read_cpgs:
            read_cpg_index +=1
            temp = []
            for pos, tag in read_cpg:
                if pos is not None and (pos > start) and (pos <= stop) and ((tag == 'Z') or (tag == 'z')):
                    # Convert from pysam's 0-based coordinate to the 1-based
                    # convention used by MethylDackel/methylKit-style reports,
                    # so downstream matrix column labels line up directly
                    # with the coordinates you'd see in a CpG report (e.g.
                    # chr1:1229803 rather than chr1:1229802).
                    temp.append((pos + 1, tag))
                    found_cpg_count += 1
                    
            output.append(temp)
        return output


    def create_matrix(self, read_cpgs):
        """
        Converted parsed reads into a pandas dataframe.

        :param read_cpgs: read CpGs generated by self.parse_reads
        :type read_cpgs: iterable

        :return: matrix methylated (1) and unmethylated (0) states
        :rtype: pd.DataFrame

        """
        series = []
        data_index = -1
        for data in read_cpgs:
            data_index += 1
            positions = []
            statues = []
            num_positions = 0
            dup_flag = False
            positions_set = set()
            for pos, status in data:
                num_positions += 1
                if pos in positions_set:
                    dup_flag = True
                positions_set.add(pos)
                
                positions.append(pos)
                statues.append(status)
            

            if dup_flag:
                pass
            else:
                if num_positions > 0:
                    series.append(pd.Series(statues, positions))

        try:
            matrix = pd.concat(series, axis=1, ignore_index=True)
        except BaseException as e:
            raise ValueError("Empty matrix")

        matrix = matrix.replace('Z', 1)
        matrix = matrix.replace('z', 0)

        return matrix.T


    def fix_read_overlap(self, full_reads, read_cpgs):
        """Takes pysam reads and read_cpgs generated during parse reads and removes any
        overlap between read1 and read2. If possible it also stitches read1 and read2 together to create
        a super read.

        :param full_reads: set of reads generated by self.parse_reads()
        :param read_cpgs: todoo
        :return: A list in the same format as read_cpgs input, but corrected for paired read overlap
        """
        # data for return
        fixed_read_cpgs = []
        # Combine raw reads and extracted tags
        combined = []
        for read, state in zip(full_reads, read_cpgs):
            combined.append((read, state))

        # Get names of all the reads present
        # query_names = [x.query_name for x in full_reads]
        query_names = []
        for x in full_reads:
            if (not (x.query_name in self.skipped_reads)) and (self.query_count_hash[x.query_name]<=2):
                query_names.append(x.query_name)

        # Match paired reads by query_name
        tally = defaultdict(list)
        for i, item in enumerate(query_names):
            tally[item].append(i)

        for key, value in sorted(tally.items()):
            # A pair exists, process it
            if len(value) == 2:
                # Set read1 and read2 correctly
                if combined[value[0]][0].is_read1:
                    read1 = combined[value[0]]
                    read2 = combined[value[1]]

                elif combined[value[1]][0].is_read1:
                    read1 = combined[value[1]]
                    read2 = combined[value[0]]

                # both reads have same value, this shouldn't be. Drop one completely, dont
                # bother with overlap
                elif combined[value[0]][0].is_read1 == combined[value[1]][0].is_read1:
                    fixed_read_cpgs.append(combined[value[0]][1])
                    continue


                else:
                    raise AttributeError("Could not determine read 1 or read 2")

                # Find amount of overlap
                amount_overlap = 0
                r1_bps = [x[0] for x in read1[1]]
                r2_bps = [x[0] for x in read2[1]]

                if min(r1_bps) < min(r2_bps):
                    trim_direction = 5
                    for bp in r2_bps:
                        if bp and bp in r1_bps:
                            amount_overlap += 1
                else:
                    trim_direction = 3
                    for bp in r1_bps:
                        if bp and bp in r2_bps:
                            amount_overlap += 1

                # remove the overlap by trimming or discarding
                if amount_overlap == len(read2[1]):
                    # discard read 2, only append read 1
                    fixed_read_cpgs.append(read1[1])
                else:
                    # trim overlap
                    if trim_direction == 5:
                        new_read2_cpgs = read2[1][amount_overlap:]
                    elif trim_direction == 3:
                        new_read2_cpgs = read2[1][:-amount_overlap]
                    # stitch together read1 and read2
                    read1[1].extend(new_read2_cpgs)
                    fixed_read_cpgs.append(read1[1])

            elif len(value) == 1:
                # No pair, add to output
                fixed_read_cpgs.append(combined[value[0]][1])

        return fixed_read_cpgs


    @staticmethod
    def correct_cpg_positions(output: list):
        """
        For some reason, Bismark alignment produces instances where a CpG site location is incorrect by 1 bp, even
        after accounting for DNA strand alignmment. This function fixes this. If two cpgs have positions such as 4, 5
        (which is impossible because there needs to by a G between them) this function will convert all 5s to 4s. This
        only needs to be applied to matrices which are empty after dropna() is called.

        :param output: a list of lists of tuples. The output of self.parse_reads()

        :return: list of the same style, execpt the first position in the tuple will have a corrected CpG position.

        """
        # find all cpg positions
        cpg_positions = []
        for item in output:
            if item:
                for cpg in item:
                    cpg_positions.append(cpg[0])
        cpg_positions = sorted(list(set(cpg_positions)))

        # determine corrections
        corrections = {}
        for x in range(len(cpg_positions)):
            try:
                if cpg_positions[x + 1] == cpg_positions[x] + 1:
                    corrections[cpg_positions[x + 1]] = cpg_positions[x]
            except IndexError:  # end of cpg position list
                pass

        # correct items
        corrected_output = []
        for item in output:
            corrected_item = []
            if item:
                for cpg in item:
                    if cpg[0] in corrections.keys():
                        new_cpg = (corrections[cpg[0]], cpg[1])
                        corrected_item.append(new_cpg)
                    else:
                        corrected_item.append(cpg)
                corrected_output.append(corrected_item)
            else:
                corrected_output.append(item)

        return corrected_output