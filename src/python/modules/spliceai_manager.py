#!/usr/bin/env python3

"""
SpliceAI prediction wrapper
"""

from .shared import CommandLineManager, dir_name_by_date, get_upper_dir
from collections import defaultdict
from heapq import heappop, heappush
from typing import Dict, List, Optional, Tuple

import click
import os

PYTHON_DIR: str = get_upper_dir(__file__, 2)
EXEC_SCRIPT: str = os.path.abspath(PYTHON_DIR, 'predict_with_spliceai.py')
CONTIG_SIZE_SCRIPT: str = os.path.join(PYTHON_DIR, 'get_contig_sizes.py')

DEFAULT_CHUNK_SIZE: int = 6_000_000
DEFAULT_FLANK_SIZE: int = 50_000
DEFAULT_MIN_CONTIG_SIZE: int = 500

DONOR_PLUS: str = 'spliceAiDonorPlus.bw'
DONOR_MINUS: str = 'spliceAiDonorMinus'
ACC_PLUS: str = 'spliceAiAcceptorPlus'
ACC_MINUS: str = 'spliceAiAcceptorMinus'

class SpliceAiManager(CommandLineManager):
    """
    """
    __slots__ = (
        'output', 'tmp_dir',
        'twobit', 'chunk_size', 'flank_size',
        'min_contig_size', 'round_to', 'min_prob',
        'job_num',
        'bed_dir', 'job_file'

        'v'
    )

    def __init__(
        self,
        ref_2bit: click.Path,
        output: Optional[click.Path],
        chunk_size: Optional[int],
        flank_size: Optional[int],
        min_contig_size: Optional[int],
        round_to: Optional[int],
        min_prob: Optional[float],
        job_number: Optional[int],
        parallel_strategy: Optional[str],
        nextflow_exec_script: Optional[click.Path],
        max_number_of_retries: Optional[int],
        nextflow_config_file: Optional[click.Path],
        max_parallel_time: Optional[int],
        cluster_queue_name: Optional[str],
        twobittofa_binary: Optional[click.Path],
        fatotwobit_binary: Optional[click.Path],
        wigtobigwig_binary: Optional[click.Path],
        keep_temporary_files: Optional[bool],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        project_name: str = dir_name_by_date('spliceai')
        self.output: str = (
            self._abspath(output) if output is not None else 
            self._abspath(project_name)
        )
        self.tmp_dir: str = os.path.join(self.output, f'tmp_{project_name}')
        log_file: str = os.path.join(self.output, f'{project_name}txt')
        self.set_logging(log_file=log_file)

        self.twobit: click.Path = ref_2bit
        self.upd_twobit: str = os.path.join(self.tmp_dir, 'unmasked_genome.2bit')

        self.chunk_size: int = chunk_size
        self.flank_size: int = flank_size
        self.min_contig_size: int = min_contig_size

        self.round_to: int = round_to
        self.min_prob: float = min_prob

        self.job_num: int = job_number

        self.bed_dir: str = os.path.jon(self.tmp_dir, 'bed_input')
        self.job_list: str = os.path.join(self.tmp_dir, 'joblist')
        self._mkdir(self.bed_dir)
        
        self._mkdir(self.output)
        self._mkdir(self.tmp_dir)
        self.run()
        if not keep_temporary_files:
            self._rmdir(self.tmp_dir)

    def run(self) -> None:
        """Entry point"""
        self.prepare_ref_genome()
        self.schedule_jobs()
        self.run_jobs()
        self.aggregate_jobs()

    def prepare_ref_genome(self) -> None:
        """
        Replaces assembly gaps with polyA in the input genome 
        since SpliceAi cannot process symbols other than A/C/G/T

        As in-place 2bit file modification is impossible, the original file
        is first decompressed, with non-standard nucleotide symbols 
        in the sequence part then being removed

        In the same run, the code extract contig sizes from the resulting Fasta file
        """
        cmd: str = (
            'set -eu; set -o pipefail; '
            f'{self.twobittofa_binary} {self.twobit} stdout | '
            f'sed \'/^>/!s/[BD-FH-SU-Z]/A/Ig\' | tee {self.tmp_fa} |'
            f'{CONTIG_SIZE_SCRIPT} - {self.chrom_sizes}'
        )
        _ = self._exec(cmd, 'Genome unmasking & contig size retrieval failed:')
        
    def schedule_jobs(self) -> None:
        """
        Schedules SpliceAI jobs for cluster execution. 
        Individual commands are grouped in {self.n_jobs} job bins. Jobs are equilibrated 
        using longest-processing-time-first (LPT) algorithm with sequence lengths 
        as proxies for execution time per command
        """
        bin2chunks: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
        job_heap: List[Tuple[int, int]] = [(0, i) for i in range(self.job_num)]
        with open(self.chrom_sizes, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != 2:
                    self._die(
                        (
                            'Improper formatting in the chromosome size file at line %s; '
                            'expected 2 columns, got % i'
                        ) % (i, len(data))
                    )
                contig, length = data
                ## ignore contigs shorter than the minimal threshold
                if length >= self.min_contig_size:
                    continue
                ## each contig is split into chunks of self.chunk_size length each + 2 * self.flank_size
                if length >= self.chunk_size:
                    ## longer sequences are split into individual chunks
                    chunk_num: int = length // self.chunk_size
                    for s in range(0, chunk_num):
                    # for s in range(0, length, self.chunk_size):
                        start: int = s * self.chunk_size
                        flanked_start: int = max(0, start - self.flank_size)
                        end: int = min(length, start + self.chunk_size)
                        flanked_end: int = min(length, end + self.flank_size)
                        chunk_size: str = end - start
                        lightest_bin: Tuple[int, int] = heappop(job_heap)
                        total_bin_size, bin_id = lightest_bin
                        name: str = f'{contig}:{start}-{end}'
                        bin2chunks[bin_id].append((contig, flanked_start, flanked_end, name))
                        total_bin_size += chunk_size
                        heappush((total_bin_size, bin_id))
                    ## if there is a trailing terminal portion, 
                    ## flank it from the left side and push it as well
                    if self.chunk_size * chunk_num < length:
                        start: int = chunk_num * self.chunk_size
                        flanked_start: int = max(0, start - self.flank_size)
                        end: int = length
                        chunk_size = end - start
                        lightest_bin: Tuple[int, int] = heappop(job_heap)
                        total_bin_size, bin_id = lightest_bin
                        name: str = f'{contig}:{flanked_start}-{end}'
                        bin2chunks[bin_id].append((contig, start, end))
                        total_bin_size += chunk_size
                        heappush((total_bin_size, bin_id))
                else:
                    ## just push the whole chunk
                    lightest_bin: Tuple[int, int] = heappop(job_heap)
                    total_bin_size, bin_id = lightest_bin
                    name: str = f'{contig}:{0}-{length}'
                    bin2chunks[bin_id].append((contig, 0, length, name))
                    total_bin_size += chunk_size
                    heappush((total_bin_size, bin_id))
        ## remove the hefty list
        del job_heap
        ## write the resulting jobs
        chunk_num: int = 0
        job_list: List[str] = []
        for bucket, intervals in bin2chunks.items():
            bed_file: str = os.path.jon(self.bed_dir, f'batch{chunk_num}.bed')
            with open(bed_file, 'w') as h:
                for interval in intervals:
                    line: str = '\t'.join(map(str, interval))
                    h.write(line + '\n')
            cmd: str = (
                f'{EXEC_SCRIPT} {self.twobit} {bed_file} ' ##TODO: Provide modified genome!
                f'--round_to {self.round_to} --min_prob {self.min_prob} -o {self.tmp_dir}'
            )
            job_list.append(cmd)
        with open(self.job_list, 'w') as h:
            for job in job_list:
                h.write(job + '\n')

    def run_jobs(self) -> None:
        """
        Controls parallel job execution
        """
        pass
        
    def aggregate_jobs(self) -> None:
        """
        Aggregates individual job results and converts them into into bigWig format,
        one per each splice site and each strand
        """
        pass

