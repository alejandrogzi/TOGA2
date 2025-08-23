"""
Reference input annotation preparation functionality
"""

__author__ = 'Yury V. Malovichko'
__credits__ = 'Michael Hiller'
__year__ = '2025'

from .constants import Headers
from contextlib import nullcontext
from .shared import (
    CommandLineManager, dir_name_by_date, 
    get_upper_dir, hex_dir_name, parse_single_column
)
from .filter_ref_bed import (
    consistent_name, 
    CONTIG_REJ_REASON, FRAME_REJ_REASON, 
    NAME_REJ_REASON, NON_CODING_REJ_REASON 
)
from shutil import which
from typing import Dict, List, Optional, Tuple, Union

import click
import os

DEFAULT_PREFIX: str = 'TOGA2_ref_annotation'
LEVEL: str = 'level'
TRANSCRIPT: str = 'TRANSCRIPT'
ILLEGAL_NAME: str = 'ILLEGAL_NAME'
REJECTED_CONTIG: str = 'REJECTED_CONTIG'
NON_CODING: str = 'NON_CODING'
OUT_OF_FRAME: str = 'OUT_OF_FRAME'
NUMERIC_FIELDS: Tuple[int] = (1, 2, 6, 7, 9, 10, 11)
U12: str = 'U12'
CANON_SITES: str = ('GT-AG', 'GC-AG')
MIN_INTRON_LENGTH_FOR_PROFILES: int = 70 ## TODO: Check

TOGA2_ROOT: str = get_upper_dir(__file__, 4)
DEFAULT_TWOBITTOFA: str = os.path.join(TOGA2_ROOT, 'bin', 'twoBitToFa')
DEFAULT_BED2FRACTION: str = os.path.join(
    TOGA2_ROOT, 'src', 'rust', 'target', 'release', 'bed12ToFraction'
)

REJ_GENE: str = (
    'GENE\t{}\t0\tNo (valid) transcripts found in the reference annotation\tZERO_TRANSCRIPT_INPUT\tN'
)
ORPHAN_TR: str = (
    'TRANSCRIPT\t{}\t0\tNo corresponding gene found in the isoform file\tZERO_GENE_INPUT\tN'
)

class InputProducer(CommandLineManager):
    """
    Main class for input annotation preparation
    """
    __slots__ = (
        'v', 'twobit', 'annot', 'isoforms', 
        'output', 'filtered_annotation', 'filtered_isoforms',
        'rejection_log', 'tr2annot', 
        'rejected_transcripts', 'rejected_lines'
        'cesar_profiles',
        'intronic', 'ic_cores',
        'twobittofa_binary', 'bed2fraction_binary',
        'intron_file', 'all_intron_bed',
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
        disable_intron_classificaiton: Optional[bool] = False,
        disable_cesar_profiles: Optional[bool] = False,
        intronic_binary: Optional[Union[click.Path,  None]] = None,
        intronic_cores: Optional[int] = 1,
        twobittofa_binary: Optional[Union[click.Path, None]] = None,
        bed2fraction_binary: Optional[Union[click.Path, None]] = None,
        min_intron_length_cesar: Optional[int] = MIN_INTRON_LENGTH_FOR_PROFILES
    ) -> None:
        self.v = True
        self.set_logging('annotation_producer')

        self.twobit: click.Path = ref_2bit
        self.annot: click.Path = ref_annot
        self.isoforms: Union[click.Path, None] = ref_isoforms
        output: str = (
            output if output is not None else hex_dir_name(DEFAULT_PREFIX)
        )

        self._mkdir(output)
        self.filtered_annotation: str = os.path.join(output, 'toga.transcripts.bed')
        self.filtered_isoforms: str = os.path.join(output, 'toga.isoforms.tsv')
        self.rejection_log: str = os.path.join(output, 'rejected_items.tsv')

        self.rejected_transcripts: List[str] = []
        self.rejected_lines: List[str] = []

        self.check_annotation(contigs, excluded_contigs)
        if self.isoforms is not None:
            self.all_transcripts: List[str] = parse_single_column(self.annot)
            self.check_isoforms()
        self.write_annotation()
        self.write_rejection_log()
        if not disable_intron_classificaiton:
            self.intronic: Union[click.Path, None] = self.check_intronic_binary(intronic_binary)
            self.ic_cores: int = intronic_cores
            self.twobittofa_binary: str = (
                twobittofa_binary if twobittofa_binary is not None else DEFAULT_TWOBITTOFA
            )
            self.bed2fraction_binary: str = (
                bed2fraction_binary if bed2fraction_binary is not None else DEFAULT_BED2FRACTION
            )
            self.intron_file: str = os.path.join(output, 'toga.U12introns.bed')
            self.all_intron_bed: str = os.path.join(output, 'all_introns.bed')
            self.intron_classifier()
        if not (disable_intron_classificaiton or disable_cesar_profiles):
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
                name: str = data[4]
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
                if frame_length % 3 and not self.no_frame_filter:
                    out_of_frame.append(name)
                    continue
                self.tr2annot[name] = line
        if illegal_name:
            self._echo(
                (
                    'The following transcripts were filtered out '
                    'due to illegal symbols used in their names:\n\t%s'
                ) % ','.join(illegal_name),
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
                ) % ','.join(rejected_contigs),
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
                ) % ','.join(non_coding),
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
                    'due to shifted reading frame'
                ) % ','.join(out_of_frame),
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
        trs_found: List[str] = []
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
                trs_found.append(tr)
        rejected_genes: List[str] = []

        ## write the remaining isoforms to the output file
        with open(self.filtered_isoforms, 'w') as h:
            for gene, trs in gene2trs:
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
                ) % ','.join(rejected_genes),
                'warning'
            )
            self.rejected_lines.extend(
                [REJ_GENE.format(x) for x in rejected_genes]
            )
        ## report transcripts for which genes were not found in the isoform file
        rejected_transcripts: List[str] = [
            x for x in self.all_transcripts if x not in trs_found
        ]
        if rejected_transcripts:
            self._to_log(
                (
                    'The following transcripts were removed from the '
                    'annotation because they were not mapped to any gene '
                    'in the isoform file'
                ) % ','.join(rejected_transcripts),
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
        with open(self.filtered_annotation, 'w') as h:
            for line in self.tr2annot.values():
                h.write(line + '\n')

    def write_rejection_log(self) -> None:
        """Writes the rejected items to the file"""
        if not self.rejected_lines:
            return
        with open(self.rejection_log, 'w') as h:
            h.write(Headers.REJ_LOG_HEADER + '\n')
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
        ## create temporary directory for intronIC meta
        tmp_dir: str = os.path.join('output', dir_name_by_date('_tmp_intronIC'))
        self._mkdir(tmp_dir)
        ## temporarily decompress the reference .2bit genome file
        self._to_log('Converting .2bit genome into Fasta format')
        genome_fasta: str = os.path.join(tmp_dir, 'genome.fa')
        decompr_cmd: str = f'{self.twobittofa_binary} {self.twobit} {genome_fasta}'
        _ = self._exec(decompr_cmd, 'twoBitToFa conversion failed:')
        ## get the intron Bed6 file
        self._to_log('Extracting the unique coding introns')
        raw_intron_file: str = os.path.join(tmp_dir, 'raw_introns.bed')
        intron_bed_cmd: str = (
            f'{self.bed2fraction_binary} -i {self.annot} -o {raw_intron_file} '
            '-m cds -b'
        )
        _ = self._exec(intron_bed_cmd, 'Intron Bed6 extraction failed:')
        ## select unique introns and anonimize them
        ## TODO: Must be merged with the previous command after writing to stdout is fixed in bed12ToFraction
        intron_bed: str = os.path.join(tmp_dir, 'intronic_input.bed')
        uniq_cmd: str = (
            f'cut -f1-3,5,6 {raw_intron_file} | sort -u | '
            'awk \'BEGIN{OFS="\t"}{print $1,$2,$3,"intron"NR,$4,$5}\' > ',
            intron_bed
        )
        _ = self._exec(uniq_cmd, 'Intron deduplication failed:')
        ## run intronIC
        self._to_log('Running intronIC')
        ic_output: str = os.path.join(tmp_dir, 'output')
        intronic_cmd: str = (
            f'{self.intronic} -g {genome_fasta} -b {intron_bed} -n {ic_output} '
            f'-p {self.ic_cores} --no_nc_ss_adjustment'
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
                intron2coords[name] = f'{data[0]}\t{data[1]}\t{data[2]}\t{{}}\t{score}\t{data[5]}'
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
                name: str = data[0].split(';')
                dinuc: str = data[2]
                intron_class: str = data[12]
                upd_name: str = f'{name}_{intron_class}_{dinuc}'
                out_line: str = intron2coords[name].format(upd_name)
                if intron_class == U12 or dinuc in CANON_SITES:
                    oh.write(out_line + '\n')
                if not self.disable_cesar_profiles:
                    ah.write(out_line + '\n')


    def cesar_profiles(self) -> None:
        """
        Given the intronIC classification for filtered annotation reference introns
        """