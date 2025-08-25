"""
Reference input annotation preparation functionality
"""

__author__ = 'Yury V. Malovichko'
__credits__ = 'Michael Hiller'
__year__ = '2025'

from collections import defaultdict
from .constants import Headers
from contextlib import nullcontext
from .shared import (
    CommandLineManager, dir_name_by_date, 
    get_upper_dir, hex_dir_name
)
from .filter_ref_bed import (
    consistent_name, 
    CONTIG_REJ_REASON, FRAME_REJ_REASON, 
    NAME_REJ_REASON, NON_CODING_REJ_REASON 
)
from shutil import which
from typing import Dict, List, Optional, Set, Tuple, Union

import click
import logging
import os

logging.basicConfig(level=logging.INFO)

DEFAULT_PREFIX: str = 'TOGA2_ref_annotation'
LEVEL: str = 'level'
TRANSCRIPT: str = 'TRANSCRIPT'
ILLEGAL_NAME: str = 'ILLEGAL_NAME'
REJECTED_CONTIG: str = 'REJECTED_CONTIG'
NON_CODING: str = 'NON_CODING'
OUT_OF_FRAME: str = 'OUT_OF_FRAME'
NUMERIC_FIELDS: Tuple[int] = (1, 2, 6, 7, 9, 10, 11)
U12: str = 'U12'
U2: str = 'U2'
CANON_SITES: Tuple[str] = ('GT-AG', 'GC-AG')
U12_CANON_SITES: str = 'GT-AG'
MIN_INTRON_LENGTH_FOR_CLASSIFICATION: int = 30
MIN_INTRON_LENGTH_FOR_PROFILES: int = 70
ENTRY_START: str = '>'
ACC: str = 'acceptor'
DONOR: str = 'donor'
ACC_PROFILE_LEN: int = 22
DONOR_PROFILE_LEN: int = 6
ACC_START_POS: int = -22
DONOR_START_POS: int = 1
NUCS: Tuple[str] = ('A', 'T', 'C', 'G')
N: str = 'N'
CANON: str = 'canon'
NONCANON: str = 'nonCanon'

TOGA2_ROOT: str = get_upper_dir(__file__, 4)
DEFAULT_TWOBITTOFA: str = os.path.join(TOGA2_ROOT, 'bin', 'twoBitToFa')
DEFAULT_BED2FRACTION: str = os.path.join(
    TOGA2_ROOT, 'src', 'rust', 'target', 'release', 'bed12ToFraction'
)
PROFILE_DIR: str = 'CESAR2.0_profiles'
EQUI_ACC: str = 'equiprobable_acceptor.tsv'
EQUI_DONOR: str = 'equiprobable_donor.tsv'

REJ_GENE: str = (
    'GENE\t{}\t0\tNo (valid) transcripts found in the reference annotation\tZERO_TRANSCRIPT_INPUT\tN'
)
ORPHAN_TR: str = (
    'TRANSCRIPT\t{}\t0\tNo corresponding gene found in the isoform file\tZERO_GENE_INPUT\tN'
)
EXTRACTION_ERR_MSG: str = 'ERROR: twoBitToFa call failed'

class InputProducer(CommandLineManager):
    """
    Main class for input annotation preparation
    """
    __slots__ = (
        'v', 'twobit', 'annot', 'isoforms', 
        'disable_intron_classification', 'disable_cesar_profiles',
        'output', 'filtered_annotation', 'filtered_isoforms',
        'rejection_log', 'tr2annot', 
        'rejected_transcripts', 'rejected_lines',
        'intronic', 'ic_cores',
        'twobittofa_binary', 'bed2fraction_binary',
        'tmp_dir', 'intron_file', 'all_intron_bed',
        'min_intron_length_intronic', 'min_intron_length_cesar',
        'intron2class', 'intron2coords', 'profiles',
        'profile_dir',
        'keep_tmp'
    )

    def __init__(
        self, 
        ref_2bit: click.Path,
        ref_annot: click.Path,
        ref_isoforms: Optional[click.Path] = None,
        output: Optional[Union[click.Path, None]] = None,
        disable_transcript_filtering: Optional[bool] = False,
        contigs: Optional[Union[str, None]] = None,
        excluded_contigs: Optional[Union[str, None]] = None,
        disable_intron_classification: Optional[bool] = False,
        disable_cesar_profiles: Optional[bool] = False,
        intronic_binary: Optional[Union[click.Path,  None]] = None,
        intronic_cores: Optional[int] = 1,
        min_intron_length_intronic: Optional[int] = MIN_INTRON_LENGTH_FOR_CLASSIFICATION,
        twobittofa_binary: Optional[Union[click.Path, None]] = None,
        min_intron_length_cesar: Optional[int] = MIN_INTRON_LENGTH_FOR_PROFILES,
        keep_temporary: Optional[bool] = False
    ) -> None:
        self.v: bool = True
        self.set_logging()

        self.twobit: click.Path = ref_2bit
        self.annot: click.Path = ref_annot
        self.isoforms: Union[click.Path, None] = ref_isoforms
        self.output: str = (
            output if output is not None else hex_dir_name(DEFAULT_PREFIX)
        )
        self.keep_tmp: bool = keep_temporary

        self._mkdir(self.output)
        self.filtered_annotation: str = os.path.join(self.output, 'toga.transcripts.bed')
        self.filtered_isoforms: str = os.path.join(self.output, 'toga.isoforms.tsv')
        self.rejection_log: str = os.path.join(self.output, 'rejected_items.tsv')

        self.tr2annot: Dict[str, str] = {}
        self.rejected_transcripts: List[str] = []
        self.rejected_lines: List[str] = []

        ## step 1: annotation file check
        self.check_annotation(contigs, excluded_contigs)

        ## step 2, optional: isoform file check, potential further annotation filtering 
        if self.isoforms is not None:
            self.check_isoforms()

        ## write the results for steps 1 and 2
        self.write_annotation()
        self.write_rejection_log()
        self.disable_intron_classification: bool = disable_intron_classification
        self.disable_cesar_profiles: bool = disable_cesar_profiles

        ## step 3: intron classification, U12 input file preparation
        if not disable_intron_classification:
            self.tmp_dir: str = os.path.join(self.output, dir_name_by_date('tmp.intronIC'))
            self.intronic: Union[click.Path, None] = self.check_intronic_binary(intronic_binary)
            self.ic_cores: int = intronic_cores
            self.twobittofa_binary: str = (
                twobittofa_binary if twobittofa_binary is not None else DEFAULT_TWOBITTOFA
            )
            self.bed2fraction_binary: str = DEFAULT_BED2FRACTION
            self.min_intron_length_intronic: int = min_intron_length_intronic
            if not disable_cesar_profiles:
                self.intron2class: Dict[str, Tuple[str, str]] = {}
                self.intron2coords: Dict[str, Tuple[str, int, int, bool]] = {}
                self.min_intron_length_cesar: int = min_intron_length_cesar
            self.intron_file: str = os.path.join(self.output, 'toga.U12introns.bed')
            self.all_intron_bed: str = os.path.join(self.tmp_dir, 'all_introns.bed')
            self._mkdir(self.tmp_dir)
            self.intron_classifier()

        ## step 4: CESAR2 profile generation; will not shoot if step 3 is disabled
        if not (disable_intron_classification or disable_cesar_profiles):
            self.profiles: Dict[Tuple[str, bool, str], Dict[int, Dict[str, int]]] = (
                defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
            )
            self.profile_dir: str = os.path.join(self.output, PROFILE_DIR)
            self._mkdir(self.profile_dir)
            self.generate_cesar_profiles()

    def check_annotation(
        self, contigs: Union[str, None], excluded_contigs: Union[str, None]
    ) -> None:
        """
        Filters reference annotation by the following criteria:
        * All transcripts in the final annotation must be 
        """
        ## TODO: Ideally copy the code here and modify as needed;
        ## a bit of silly code repetition, but at least no need to parse the rejection log
        illegal_name: List[str] = []
        rejected_contigs: List[str] = []
        non_coding: List[str] = []
        out_of_frame: List[str] = []
        with open(self.annot, 'r') as h:
            for i, line in enumerate(h, start=1):
                line = line.strip()
                data: List[str] = line.split('\t')
                if not data or not data[0]:
                    continue
                if len(data) != 12 :
                    self._die(
                        (
                            'Improper formatting at reference annotation file line %i; '
                            'expected 12 fields, got %i'
                        ) % (i, len(data))
                    )
                if any(x.strip() == '' for x in data):
                    self._die(
                        (
                        'Improper formatting at reference annotation file line %i; '
                        'empty fields encountered'
                        ) % i
                    )
                for field in NUMERIC_FIELDS:
                    if not data[field].replace(',', '').isdigit():
                        self._die(
                            (
                            'Improper formatting at reference annotation file line %i; '
                            'field %i contains non-numeric data'
                            ) % (i, field)
                        )
                name: str = data[3]
                ## remove the transcripts with improperly formatted names
                if not consistent_name(name):
                    illegal_name.append(name)
                    continue
                chrom: str = data[0]
                ## if entries were restricted to specific contigs,
                ## apply the respective filters
                if contigs and chrom not in contigs:
                    rejected_contigs.append(name)
                    continue
                if excluded_contigs and chrom in excluded_contigs:
                    rejected_contigs.append(name)
                    continue
                ## check coding sequence presence and frame intactness
                thin_start: int = int(data[1])
                # thin_end: int = int(data[2])
                cds_start: int = int(data[6])
                cds_end: int = int(data[7])
                if cds_end < cds_start:
                    self._die(
                        (
                            'Improper formatting at reference annotation file line %i; '
                            'coding sequence start coordinate greated than the start coordinate'
                        )
                    )
                if cds_start == cds_end:
                    non_coding.append(name)
                    continue
                ## iterate over exon entries to infer the CDS length
                frame_length: int = 0
                sizes: List[int] = list(map(int, data[10].split(',')[:-1]))
                starts: List[int] = list(map(int, data[11].split(',')[:-1]))
                for start, size in zip(starts, sizes):
                    start += thin_start
                    end: int = start + size
                    if start < cds_start:
                        if end > cds_start:
                            size -= (cds_start - start)
                        else:
                            continue
                    if end > cds_end:
                        if start < cds_end:
                            size -= (end - cds_end)
                        else:
                            continue
                    frame_length += size
                if frame_length % 3:# and not self.no_frame_filter:
                    out_of_frame.append(name)
                    continue
                self.tr2annot[name] = line
        if illegal_name:
            self._to_log(
                (
                    'The following transcripts were filtered out '
                    'due to illegal symbols used in their names:\n\t%s'
                ) % '\n\t'.join(illegal_name),
                'warning'
            )
            self.rejected_transcripts.extend(illegal_name)
            self.rejected_lines.extend(
                [NAME_REJ_REASON.format(x) for x in illegal_name]
            )
        if rejected_contigs:
            self._to_log(
                (
                    'The following transcripts were filtered out '
                    'due to their location in deprecated contigs:\n\t%s'
                ) % '\n\t'.join(rejected_contigs),
                'warning'
            )
            self.rejected_transcripts.extend(rejected_contigs)
            self.rejected_lines.extend(
                [CONTIG_REJ_REASON.format(x) for x in rejected_contigs]
            )
        if non_coding:
            self._to_log(
                (
                    'The following transcripts were filtered out '
                    'due to absence of coding sequence:\n\t%s'
                ) % '\n\t'.join(non_coding),
                'warning'
            )
            self.rejected_transcripts.extend(non_coding)
            self.rejected_lines.extend(
                [NON_CODING_REJ_REASON.format(x) for x in non_coding]
            )
        if out_of_frame:
            self._to_log(
                (
                    'The following transcripts were filtered out '
                    'due to shifted reading frame:\n\t%s'
                ) % '\n\t'.join(out_of_frame),
                'warning'
            )
            self.rejected_transcripts.extend(out_of_frame)
            self.rejected_lines.extend(
                [FRAME_REJ_REASON.format(x) for x in out_of_frame]
            )
        ## proceed further

    def check_isoforms(self) -> None:
        """
        Filters the reference isoform (gene-to-transcript) mapping file
        by removing genes whose transcripts (isoforms) were discarded 
        at the annotation filter step. 
        In parallel, all the transcripts recorded at the annotation filter step 
        and missing gene mapping in the isoforms file are further excluded 
        from the final annotation. 
        """
        gene2trs: Dict[str, List[str]] = {}
        trs_found: Set[str] = set()
        with open(self.isoforms, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.rstrip().split('\t')
                if not data or not data[0]:
                    continue
                if len(data) != 2:
                    self._die(
                        (
                            'Improper formatting at isoforms file line %i; '
                            'expecting 2 columns, got %i'
                        ) % (i, len(data))
                    )
                gene, tr = data
                if gene not in gene2trs:
                    gene2trs[gene] = []
                if tr in self.rejected_transcripts:
                    continue
                gene2trs[gene].append(tr)
                trs_found.add(tr)
        rejected_genes: List[str] = []

        ## write the remaining isoforms to the output file
        with open(self.filtered_isoforms, 'w') as h:
            for gene, trs in gene2trs.items():
                if not trs:
                    rejected_genes.append(gene)
                    continue
                for tr in trs:
                    h.write(f'{gene}\t{tr}\n')

        ## report genes which ended up having no transcripts
        if rejected_genes:
            self._to_log(
                (
                    'The following genes were removed from the '
                    'isoforms file because all the respective '
                    'transcripts were removed from the annotation:\n\t%s'
                ) % '\n\t'.join(rejected_genes),
                'warning'
            )
            self.rejected_lines.extend(
                [REJ_GENE.format(x) for x in rejected_genes]
            )
        ## report transcripts for which genes were not found in the isoform file
        rejected_transcripts: List[str] = [
            x for x in self.tr2annot if x not in trs_found
        ]
        if rejected_transcripts:
            self._to_log(
                (
                    'The following transcripts were removed from the '
                    'annotation because they were not mapped to any gene '
                    'in the isoform file'
                ) % '\n\t'.join(rejected_transcripts),
                'warning'
            )
            self.rejected_lines.extend(
                [ORPHAN_TR.format(x) for x in rejected_transcripts]
            )
            self.tr2annot = {
                k:v for k,v in self.tr2annot.items() if k not in rejected_transcripts
            }

    def write_annotation(self) -> None:
        """Writes the filtered reference annotation to the file"""
        if not self.tr2annot:
            self._die(
                'All transcripts were filtered out for various reasons'
            )
        self._to_log('Writing the filtered annotation to %s' % self.filtered_annotation)
        with open(self.filtered_annotation, 'w') as h:
            for line in self.tr2annot.values():
                h.write(line + '\n')

    def write_rejection_log(self) -> None:
        """Writes the rejected items to the file"""
        if not self.rejected_lines:
            return
        self._to_log('Writing the rejected items to %s' % self.rejection_log)
        with open(self.rejection_log, 'w') as h:
            h.write(Headers.REJ_LOG_HEADER)
            for line in self.rejected_lines:
                h.write(line + '\n')

    def check_intronic_binary(self, binary: Union[str, None]) -> str:
        """
        Checks intronIC binary availability and execution permissions. 

        Since the input is provided from the click Python CLI, the existence 
        of non-empty option value is basically guaranteed at this point; 
        nevertheless, the code further checks whether the provided intronIC 
        instance is executable.

        If no value is provided, the method searches for intronIC availability in $PATH. 
        Once it is found, the execution permissions are further ensured.

        The method throws error if no (executable) intronIC instance is found.
        """
        if binary is not None:
            self._to_log('Testing intronIC binary at %s' % binary)
            if os.access(binary, os.X_OK):
                self._to_log(
                    'The provided binary is executable; using the stated intronIC instance'
                )
                return binary
            else:
                self._to_log(
                    (
                        'intronIC binary at %s does not seem executable; '
                        'looking for alternatives in $PATH'
                    ) % binary,
                    'warning'
                )
        else:
            self._to_log(
                (
                    'No intronIC executable was provided; '
                    'looking for alternatives in $PATH'
                )
            )
        intronic_in_path: Union[str, None] = which('intronIC')
        if intronic_in_path is not None:
            self._to_log(
                'Found intronIC instance at %s; checking the execution permissions' % intronic_in_path
            )
            if os.access(binary, os.X_OK):
                self._to_log(
                    'The found binary is executable; using the stated intronIC instance'
                )
                return intronic_in_path
            self._die(
                (
                    'The intronIC binary found in $PATH at %s is not executable; '
                    'check your $PATH or provide a valid intronIC instance'
                ) % intronic_in_path
            )
        self._die(
            (
                'No intronIC binary found in $PATH; '
                'check your $PATH or provide a valid intronIC instance'
            )
        )


    def intron_classifier(self) -> None:
        """
        Classifies the introns in the filtered annotation, separating them
        into U2 and U12 classes according to the intronIC predictions and adding 
        terminal dinucleotide data for further classification into canonical (GT-AG)
        and non-canonical (any other dinucleotide combination) subclasses.

        The resulting file is a Bed6 file of class- and dinucleotide-annotated reference introns. 
        This file can be further used with TOGA2 runs provided with the --ut12_file/-u2 option. 
        If CESAR2 profile annotation is requested, the same predictions are further used 
        to generate reference-specific profiles. 
        """
        ## temporarily decompress the reference .2bit genome file
        self._to_log('Converting .2bit genome into Fasta format')
        genome_fasta: str = os.path.join(self.tmp_dir, 'genome.fa')
        decompr_cmd: str = f'{self.twobittofa_binary} {self.twobit} {genome_fasta}'
        _ = self._exec(decompr_cmd, 'twoBitToFa conversion failed:')
        ## get the intron Bed6 file
        self._to_log('Extracting the unique coding introns')
        raw_intron_file: str = os.path.join(self.tmp_dir, 'raw_introns.bed')
        intron_bed_cmd: str = (
            f'{self.bed2fraction_binary} -i {self.annot} -o {raw_intron_file} '
            '-m cds -n -b'
        )
        _ = self._exec(intron_bed_cmd, 'Intron Bed6 extraction failed:')
        ## select unique introns and anonimize them
        ## TODO: Must be merged with the previous command after writing to stdout is fixed in bed12ToFraction
        intron_bed: str = os.path.join(self.tmp_dir, 'intronic_input.bed')
        uniq_cmd: str = (
            f'cut -f1-3,5,6 {raw_intron_file} | sort -u | '
            'awk \'BEGIN{OFS="\\t"}{print $1,$2,$3,"intron"NR,0,$5}\' > '
            f'{intron_bed}'
        )
        _ = self._exec(uniq_cmd, 'Intron deduplication failed:')
        ## run intronIC
        self._to_log('Running intronIC')
        ic_output: str = os.path.join(self.tmp_dir, 'output')
        intronic_cmd: str = (
            f'{self.intronic} -g {genome_fasta} -b {intron_bed} -n {ic_output} '
            f'-p {self.ic_cores}  --min_intron_len {self.min_intron_length_intronic} '
            '--no_nc_ss_adjustment --no_abbreviate'
        )
        _ = self._exec(intronic_cmd, 'intronIC run failed:')
        ## now, prepare the final file
        ## parse the output bed file
        ## score field contains probability, in per cent value
        ## for the final Bed file, those values are multiplied by 10 (to confine within 0<=x<=1000)
        ## and rounded to the closest integer
        out_bed_file: str = f'{ic_output}.bed.iic'
        intron2coords: Dict[str, str] = {}
        with open(out_bed_file, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != 6:
                    self._die(
                        (
                            'Improper intronIC bed file formatting at line % i; '
                            'expected 6 fields, got %i'
                        ) % (i, len(data))
                    )
                name: str = data[3].split(';')[0]
                ## convert the probability into a Bed-compatible score
                score: int = 0 if data[4] == '.' else int(float(data[4]) * 10)
                # intron2coords[name] = f'{data[0]}\t{data[1]}\t{data[2]}\t{{}}\t{score}\t{data[5]}'
                intron2coords[name] = (data[0], data[1], data[2], '{}', str(score), data[5])
        ## all is left is intron class and terminal dinucleotides
        ## extract those from the .meta.iic file, then write the results to file
        out_meta_file: str = f'{ic_output}.meta.iic'
        with (
            open(out_meta_file, 'r') as ih, 
            open(self.intron_file, 'w') as oh,
            (nullcontext if self.disable_cesar_profiles else open(self.all_intron_bed, 'w')) as ah
        ):
            for i, line in enumerate(ih, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != 14:
                    self._die(
                        (
                            'Improper intronIC file formatting at line % i; '
                            'expected 14 fields, got %i'
                        ) % (i, len(data))
                    )
                name: str = data[0].split(';')[0]
                dinuc: str = data[2]
                intron_class: str = data[12].upper()
                upd_name: str = f'{name}_{intron_class}_{dinuc}'
                out_line: str = '\t'.join(intron2coords[name]).format(upd_name)
                ## for TOGA2 input, only U12 and non-canonical U2 introns are required
                is_canon: bool = (
                    intron_class == U2 and dinuc in CANON_SITES
                ) or (intron_class == U12 and dinuc == U12_CANON_SITES)
                if intron_class == U12 or not is_canon:
                    oh.write(out_line + '\n')
                ## if CESAR2 profile generation was requested, 
                ## save all the introns to the provisional file
                ## and store their class data
                if not self.disable_cesar_profiles:
                    ah.write(out_line + '\n')
                    chrom, start, end, _, _, strand = intron2coords[name]
                    start = int(start)
                    end = int(end)
                    strand: bool = strand == '+'
                    if end - start < self.min_intron_length_cesar:
                        continue
                    intron_key: Tuple[str, int, int, bool] = (chrom, start, end, strand)
                    num: str = name.replace('intron', '')
                    self.intron2class[num] = (intron_class, is_canon)
                    self.intron2coords[num] = intron_key

        ## if profile generation is disabled and/or not requested to leave,
        ## remove the temporary directory
        if not self.keep_tmp or self.disable_cesar_profiles:
            self._rmdir(self.tmp_dir)



    def generate_cesar_profiles(self) -> None:
        """
        Given the intronIC classification for filtered annotation reference introns
        """
        ## step 1: prepare profile sequence coordinates in Bed format
        bed_string: str = ''
        ## extract coordinates
        for num, coords in self.intron2coords.items():
            chrom, start, end, strand = coords
            ## define donor and acceptor sequence coordinates
            if strand:
                donor_start: int = start
                donor_end: int = min(donor_start + DONOR_PROFILE_LEN, end)
                acc_start: int = max(end - ACC_PROFILE_LEN, start)
                acc_end: int = end
            else:
                donor_start: int = max(end - DONOR_PROFILE_LEN, start)
                donor_end: int = end
                acc_start: int = start
                acc_end: int = min(start + ACC_PROFILE_LEN, end)
            bed_strand: str = '+' if strand else '-'
            donor_name: str = f'{num}_{DONOR}'
            acc_name: str = f'{num}_{ACC}'
            donor_bed_line: str = '\t'.join(
                map(
                    str,
                    [chrom, donor_start, donor_end, donor_name, 0, bed_strand]
                )
            )
            acc_bed_line: str = '\t'.join(
                map(
                    str,
                    [chrom, acc_start, acc_end, acc_name, 0, bed_strand]
                )
            )
            bed_string += donor_bed_line + '\n' + acc_bed_line + '\n'
        ## then extract the sequences in Fasta format
        cmd: str = f'{self.twobittofa_binary} -bed=/dev/stdin {self.twobit} stdout'
        res: str = self._exec(
                cmd,
                err_msg=EXTRACTION_ERR_MSG,
                input_=bed_string.encode('utf8')
            )
        
        ## step 2: parse the resulting Fasta and update the resulting profiles
        self._to_log('Inferring positional letter probabilities')
        header: str = ''
        seq: str = ''
        for line in res.split('\n'):
            line = line.rstrip()
            if not line:
                continue
            if line[0] == ENTRY_START:
                if header:
                    num, site_type = header.split('_')
                    seq = seq.upper().strip('\n')
                    intron_class, canon = self.intron2class[num]
                    ## update the respective CESAR2 profile
                    self._update_profile(seq, intron_class, canon, site_type, num)
                    seq = ''
                header = line[1:]
            else:
                seq += line
        if seq:
            num, site_type = header.split('_')
            seq = seq.upper()
            intron_class, canon = self.intron2class[num]
            ## again, update the respective profile for the final item
            self._update_profile(seq, intron_class, canon, site_type, num)

        ## step 3: compute the letter probabilities per position from the recorded frequencies
        ## and write the resulting profiles
        self._to_log(
            'Computing CESAR2 profiles and writing them to %s' % self.profile_dir
        )
        header: str = '\t'.join(NUCS)
        for key in self.profiles:
            intron_class, canon, site_type = key
            canon_line: str = CANON if canon else NONCANON
            filename: str = f'{canon_line}_{intron_class}_{site_type}.tsv'
            with open(os.path.join(self.profile_dir, filename), 'w') as h:
                h.write(header + '\n')
                for pos in sorted(self.profiles[key].keys()):
                    pos_sum: float = float(sum(self.profiles[key][pos].values()))
                    frequencies: List[float] = []
                    for nuc in NUCS:
                        freq: float = round(
                            self.profiles[key][pos].get(nuc, 0.0) / pos_sum,
                            3
                        )
                        frequencies.append(freq)
                    h.write('\t'.join(map(str, frequencies)) + '\n')
        
        ## generate equiprobable acceptor and donor profiles
        ## equiprobable acceptor is recommended for non-canonical U12 in mammals
        ## equiprobable donors have not been tested for any purpose yet
        ## but it's nice to have both options generated automatically
        equi_line: str = '\t'.join(map(str, [0.25] * 4)) + '\n'
        equi_acc: str = os.path.join(self.profile_dir, EQUI_ACC)
        with open(equi_acc, 'w') as h:
            h.write(header + '\n')
            for _ in range(ACC_PROFILE_LEN):
                h.write(equi_line)
        equi_donor: str = os.path.join(self.profile_dir, EQUI_DONOR)
        with open(equi_donor, 'w') as h:
            h.write(header + '\n')
            for _ in range(DONOR_PROFILE_LEN):
                h.write(equi_line)

        ## remove temporary directory unless requested to leave
        if not self.keep_tmp:
            self._rmdir(self.tmp_dir)


    def _update_profile(
        self, seq: str, intron_class: str, canon: bool, site: str, num: str
    ) -> None:
        """
        Updates the CESAR2 profile with respect to spliceosomal class, 'canonicity', and splice site
        """
        if site not in (DONOR, ACC):
            self._die('Unknown site type received as input: %s' % site)
        key: Tuple[str, bool, str] = (intron_class, canon, site)
        ## CESAR2 operates fixed profile lengths for donor and acceptor 
        exp_len: int = DONOR_PROFILE_LEN if site == DONOR else ACC_PROFILE_LEN
        ## sanity check: extracted sequence must not be longer than the profile length
        ## the reverse indicates a clear extraction error
        if len(seq) > exp_len:
            self._die(
                'Sequence %s is longer than the expected %s profile length of %i' % (
                    seq, site, exp_len
                )
            )
        ## earlier TOGA2 code threw the same error if the profile was shorter than the consensus length
        ## however, for certain clades introns <22bp are quite common
        if len(seq) < exp_len:
            self._to_log(
                'Sequence %s is shorter than the expected %s profile length of %i' % (
                    seq, site, exp_len
                ),
                'warning'
            )
        ## position numeration differs 
        start_pos: int = DONOR_START_POS if site == DONOR else ACC_START_POS
        for i in range(exp_len):
            ## WARNING: The trick was not tested in TOGA2 alpha since introns <30 bp
            ## were ignored in production runs and introns <70bp were excluded from profile preparation
            ## but the idea is to make position beyond the average intron length potentially equiprobable
            pos: int = i + start_pos
            if i >= len(seq):
                for nuc in NUCS:
                    self.profile[key][pos][nuc] + 1
                continue
            nuc: str = seq[i]
            if nuc == N:
                continue
            if nuc not in NUCS:
                self._die('Ambiguous nucleotide encountered: %s' % nuc)
            self.profiles[key][pos][nuc] += 1

    def set_logging(self) -> None:
        """Sets logging and disables logging propagation"""
        super().set_logging()
        self.logger.propagate = False

