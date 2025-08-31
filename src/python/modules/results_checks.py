"""
Sanity check and troubleshooting manager
"""

from collections import defaultdict
from dataclasses import dataclass
from .constants import Constants
from .shared import CommandLineManager, parse_single_column
from typing import (
    Dict, List, Optional, Set, Tuple, Union
)

import logging
import os
import sys


ORTH_PROB_HEADER: str = 'transcript'
ORTH_PROB_FIELDS: int = 3
MEM_FILE_HEADER: str = 'transcript'
ISOFORMS_FIELDS: int = 2
MEM_FILE_FIELDS: int = 15
REJ_LOG_FIELDS: int = 6
CHAIN_LIMIT_REASON: str = 'CHAIN_LIMIT_EXCEEDED'
LOSS_HEADER: str = 'level'
LOSS_FILE_FIELDS: int = 3
PROJECTION: str = 'PROJECTION'
TRANSCRIPT: str = 'TRANSCRIPT'
GENE: str = 'GENE'
ORTH_CLASS_HEADER: str = 't_gene'
ORTH_CLASS_FIELDS: int = 5
IGNORED_STATUSES: Tuple[str, str, str] = ('PG', 'PP', 'N')
QUERY_ISOFORMS_HEADER: str = 'query_gene'
REG: str = 'reg_'

MIN_ORTH_PERCENTAGE: float = 50.0
MAX_REJ_AT_PREPR_PERCENTAGE: float = 1 / 3
MIN_INTACT_PROJ_PERCENTAGE: float = 50.0
MIN_ONE2ONE_PERCENTAGE: float = 50.0
MAX_ONE2ZERO_PERCENTAGE: float = 25.0


BREAK_LINE: str = '#' * 100
SUMMARY_BOILERPLATE: str = """
{br}
#### TOGA2 run summary
{br}

Reference genome: {ref_2bit}
Query genome: {query_2bit}
Genome alignment chain file: {chains}
Reference annotation file: {ref_annot}
Output directory: {output}

{isoforms_line}

{br}

{orthology_class_report}

{br}

{loss_summary}

{br}

{orthology_res_report}

"""

ISOFORMS_LINE: str = 'Isoforms file: {}\n'
ISOFORMS_CAVEAT: str = """
\nNOTE: Reference isoform file was not provided. All gene level statistics are presented \
assuming that each reference transcript represents a separate gene
"""

CLASSIFICATION_BOILERPLATE: str = """
Orthology prediction statistics:
\tReference transcripts level::
\t\t#reference transcripts subjected to classification: {}
\t\t#reference transcripts with >=1 predicted ortholog: {} ({}%)
\t\t#reference transcripts with 1 predicted ortholog: {} ({}%)
\t\t#reference transcripts with no predicted orthologs: {} ({}%)
\t\t#reference transcripts with no classifiable projections: {} ({}%)
\t\t#predicted processed pseudogene/retrogene projections {}:
\tReference genes level:
\t\t#reference genes with >=1 predicted ortholog: {} ({}%)
\t\t#reference genes with 1 predicted ortholog: {} ({}%)
\t\t#reference genes with no predicted orthologs: {} ({}%)
"""
LOSS_SUMMARY_BOILERPLATE: str = """
Gene loss summary statistics:
\tloss classes considered for assessing gene presence: {}
\tQuery projections level:
\t\t#Fully Intact (FI) - {} ({}%)
\t\t#Intact (I) - {} ({}%)
\t\t#Partially Intact (PI) - {} ({}%)
\t\t#Uncertain Losses (UL) - {} ({}%)
\t\t#Lost (L) - {} ({}%)
\t\t#Missing (M) - {} ({}%)
\t\t#projections considered present: {} ({}%)
\t\t#projections considered lost/missing: {} ({}%)
\tReference transcripts level:
\t\t#Fully Intact (FI) - {} ({}%)
\t\t#Intact (I) - {} ({}%)
\t\t#Partially Intact (PI) - {} ({}%)
\t\t#Uncertain Losses (UL) - {} ({}%)
\t\t#Lost (L) - {} ({}%)
\t\t#Missing (M) - {} ({}%)
\t\t#trancripts considered present: {} ({}%)
\t\t#transcripts considered lost/missing: {} ({}%)
\tReference genes level:
\t\t#Fully Intact (FI) - {} ({}%)
\t\t#Intact (I) - {} ({}%)
\t\t#Partially Intact (PI) - {} ({}%)
\t\t#Uncertain Losses (UL) - {} ({}%)
\t\t#Lost (L) - {} ({}%)
\t\t#Missing (M) - {} ({}%)
\t\t#reference genes considered present: {} ({}%)
\t\t#reference genes considered lost/missing: {} ({}%)
"""

ORTHOLOGY_BOILERPLATE: str = """
Orthology resolution:
\t#reference genes: {}
\t#query genes: {}
\t\t#with defined orthology: {} ({}%)
\t\t#lost, missing, or lacking defined orthology: {} ({}%) 
\tReference gene orthology class composition:
\t\tone2one: {} ({}%)
\t\tone2many: {} ({}%)
\t\tmany2one: {} ({}%)
\t\tmany2many: {} ({}%)
\t\tone2zero: {} ({}%)
"""

CLASSIFICATION_STUB: str = """\
{num} reference {lvl} ({fraction}%) do not have any predicted orthologous projections\
"""
CLASSIFICATION_OK: str = 'Sanity check passed: ' + CLASSIFICATION_STUB
CLASSIFICATION_WARNING: str = (
    'Sanity check failed: ' + CLASSIFICATION_STUB + 
    ', indicating potential problems with the current TOGA2 run'
)
PREPROCESSING_STUB: str = """\
Of {} projections having reached the preprocessing step, {} projections ({}%) were discarded\
"""
PREPROCESSING_OK: str = 'Sanity check passed: ' + PREPROCESSING_STUB
PREPROCESSING_WARNING: str = (
    'Sanity check failed: ' + PREPROCESSING_STUB + 
    ', indicating potential problems with the current TOGA2 run'
)
LOSS_STUB: str = """\
Of {total} orthologous {lvl} having reached the gene loss summary step, \
{num_valid} {lvl} ({fraction}%) do not belong to loss classes \
designating gene presence ({classes})\
"""
LOSS_OK: str = 'Sanity check passed: ' + LOSS_STUB
LOSS_WARNING: str = (
    'Sanity check failed: ' + LOSS_STUB + 
    ', indicating potential problems with the current TOGA2 run'
)
RESOLUTION_STUB: str = """\
Of {total} reference genes, {num} ({fraction}%) genes have {type} orthology status\
"""
RESOLUTION_OK: str = 'Sanity check passed: ' + RESOLUTION_STUB
RESOLUTION_WARNING: str = (
    'Sanity check failed: ' + RESOLUTION_STUB + 
    ', indicating potential problems with the current TOGA2 run'
)

@dataclass
class ElementaryCheckResult:
    is_warning: bool
    log_message: str
    warning_message: str


@dataclass
class SanityCheckResult:
    step: str
    messages: List[ElementaryCheckResult]


def to_perc(numerator: int, denominator: int) -> float:
    """
    Divides numerator by denominator, multiplies by 100, 
    and rounds to third digit after the dot
    """
    return round(numerator / denominator * 100, 3)


class ResultChecker(CommandLineManager):
    """
    Checks output at selected TOGA2 pipeline steps; 
    at the end of the pipeline, summarizes results and 
    """
    __slots__ = (
        'logger', 'orth_probs', 
        'feat_rej_log', 'class_rej_log',
        'memory_report', 'prepr_rej_log', 
        'loss_summary', 'orth_classes',
        'query_genes_file',
        'isoforms', 'paralogs', 'ppgenes',
        'orth_thresh', 'intact_classes',
        'all_initial_transcripts',
        'tr2orths', 'gene2tr',
        'non_orthologous', 
        'num_orth', 'fr_orth', 
        'num_o2o', 'fr_o2o', 
        'num_no_orth', 'fr_no_orth', 
        'num_no_chains', 'fr_no_chains',
        'num_ppgene', 
        'gene_num_orth', 'gene_fr_orth',
        'gene_num_o2o', 'gene_fr_o2o', 
        'gene_num_no_orth', 'gene_fr_no_orth',
        'proj2loss', 'tr2loss',
        'gene2loss', 'orth_status_counter', 
        'gene2orth_class'
    )

    def __init__(
        self,
        logger: logging.Logger,
        output: str,
        ref_isoforms: Optional[Union[str, None]] = None,
        orthology_threshold: Optional[float] = Constants.DEFAULT_ORTH_THRESHOLD,
        intact_classes: Optional[Union[List, str]] = Constants.DEFAULT_LOSS_SYMBOLS
    ) -> None:
        pass
        ## TODO: If pipeline was resumed from a certain step, 
        ## do not forget to check the results 
        self.logger: logging.Logger = logger
        self.orth_probs: str = os.path.join(
            output, 'orthology_scores.tsv'
        )
        self.feat_rej_log: str = os.path.join(
            output, 'tmp', 'rejection_logs', 'rejected_at_feature_extraction.tsv'
        )
        self.class_rej_log: str = os.path.join(
            output, 'tmp', 'rejection_logs', 'rejected_at_classification.tsv'
        )
        self.memory_report: str = os.path.join(
            output, 'meta', 'memory_requirements.tsv'
        )
        self.prepr_rej_log: str = os.path.join(
            output, 'tmp', 'rejection_logs', 'rejected_at_preprocessing.tsv'
        )
        self.loss_summary: str = os.path.join(
            output, 'loss_summary.tsv'
        )
        self.orth_classes: str = os.path.join(
            output, 'tmp', 'annotation_raw', 'orthology_classification.tsv'
        )
        self.query_genes_file: str = os.path.join(
            output, 'query_genes.tsv'
        )
        self.isoforms: Union[str, None] = ref_isoforms
        self.paralogs: str = os.path.join(
            output, 'meta', 'paralogous_projections_to_align.tsv'
        )
        self.ppgenes: str = os.path.join(
            output, 'meta', 'processed_pseudogene_projections_to_align.tsv'
        )
        self.orth_thresh: float = orthology_threshold
        if isinstance(intact_classes, str):
            if intact_classes == 'ALL':
                self.intact_classes: List[str] = Constants.ALL_LOSS_SYMBOLS
            else:
                self.intact_classes: List[str] = [
                    x.strip() for x in intact_classes.split(',') if x
                ]
        else:
            self.intact_classes: List[str] = intact_classes

        self.all_initial_transcripts: Set[str] = set()
        self.tr2orths: Dict[str, List[str]] = defaultdict(list)
        self.gene2tr: Dict[str, List[str]] = defaultdict(list)

        self.non_orthologous: Set[str] = set()

        self.num_orth: int = 0
        self.fr_orth: float = 0.0
        self.num_o2o: int = 0
        self.fr_o2o: float = 0.0
        self.num_no_orth: int = 0
        self.fr_no_orth: float = 0.0
        self.num_no_chains: int = 0
        self.fr_no_chains: float = 0.0
        self.num_ppgene: int = 0
        self.gene_num_orth: int = 0
        self.gene_fr_orth: float = 0.0
        self.gene_num_o2o: int = 0
        self.gene_fr_o2o: float = 0.0
        self.gene_num_no_orth: int = 0
        self.gene_fr_no_orth: float = 0.0

        self.proj2loss: Dict[str, int] = defaultdict(int)
        self.tr2loss: Dict[str, int] = defaultdict(int)
        self.gene2loss: Dict[str, int] = defaultdict(int)

        self.orth_status_counter: Dict[str, int] = defaultdict(int)

        self.gene2orth_class: Dict[str, int] = defaultdict(int)

    def _throw_warning(
        self, 
        text: str, 
        to_email: Optional[bool] = False
    ) -> Union[str, None]:
        """Throw warning to the provided Logger channel"""
        pass
        if to_email:
            pass ## return warning text

    def _parse_isoform_file(self) -> None:
        """Retrieves gene-to-transcripts mapping for reference genome annotation"""
        with open(self.isoforms, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().spit('\t')
                if not data:
                    continue
                if len(data) != ISOFORMS_FIELDS:
                    self.logger.critical(
                        (
                            'Improper formatting at reference isoform file line %i; '
                            'expected %i columns, got %i'
                        ) % (i, ISOFORMS_FIELDS, len(data))
                    )
                    sys.exit(1)
                gene, tr = data
                self.gene2tr[gene].append(tr)

    def _parse_non_orthologous(self) -> None:
        """Parses paralogous and retrocopy projection files"""
        if os.path.exists(self.paralogs):
            self.non_orthologous = self.non_orthologous.union(
                parse_single_column(self.paralogs)
            )
        if os.path.exists(self.ppgenes):
            self.non_orthologous = self.non_orthologous.union(
                parse_single_column(self.ppgenes)
            )

    def _get_query_gene_num(self) -> Tuple[int, int]:
        """Counts the number of query genes"""
        all_query_genes: Set[str] = set()
        if not os.path.exists(self.query_genes_file):
            self.logger.critical(
                'Reference genes file %s does not exist' % self.query_genes_file
            )
            sys.exit(1)
        with open(self.query_genes_file, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != ISOFORMS_FIELDS:
                    self.logger.critical(
                        (
                            'Improper formatting at query isoform file line %i; '
                            'expected %i columns, got %i'
                        ) % (i, ISOFORMS_FIELDS, len(data))
                    )
                    sys.exit(1)
                if data[0] == QUERY_ISOFORMS_HEADER:
                    continue
                all_query_genes.add(data[0])

        orphan: int = sum(x.startswith(REG) for x in all_query_genes)
        orth: int = len(all_query_genes) - orphan

        return (orth, orphan)

    def check_classification(self) -> SanityCheckResult:
        """
        Checks orthology prediction results. Throws a warning if: 
        a) <50% of the transcripts do not have at least one ortholog;
        b) <50% reference genes do not have any isoform with at least one ortholog
        """
        if not os.path.exists(self.orth_probs):
            self.logger.critical(
                'Orthology score fiel %i does not exist' % self.orth_probs
            )
            sys.exit(1)
        messages: List[ElementaryCheckResult] = []
        with open(self.orth_probs, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != ORTH_PROB_FIELDS:
                    self.logger.critical(
                        (
                            'Improper formatting at orthology probability file line %i; '
                            'expected %i columns, got %i'
                        ) % (i, ORTH_PROB_FIELDS, len(data)),
                    )
                    sys.exit(1)
                if data[0] == ORTH_PROB_HEADER:
                    continue
                tr: str = data[0]
                self.all_initial_transcripts.add(tr)
                chain: str = data[1]
                prob: float = float(data[2])
                proj: str = f'{tr}#{chain}'
                if prob >= self.orth_thresh:
                    self.tr2orths[tr].append(proj)
                if prob == Constants.PPGENE_PROB:
                    self.num_ppgene += 1

        discarded_at_classification: Set[str] = set()
        if os.path.exists(self.feat_rej_log):
            with open(self.feat_rej_log, 'r') as h:
                for i, line in enumerate(h, start=1):
                    data: List[str] = line.strip().split('\t')
                    if not data:
                        continue
                    if len(data) != REJ_LOG_FIELDS:
                        self.logger.critical(
                            (
                                'Improper formatting at feature extraction rejection log file line %i; '
                                'expected %i columns, got %i'
                            ) % (i, REJ_LOG_FIELDS, len(data)),
                        )
                        sys.exit(1)
                    if data[0] != TRANSCRIPT:
                        continue
                    discarded_at_classification.add(data[1])
        if os.path.exists(self.class_rej_log):
            with open(self.class_rej_log, 'r') as h:
                for i, line in enumerate(h, start=1):
                    data: List[str] = line.strip().split('\t')
                    if not data:
                        continue
                    if len(data) != REJ_LOG_FIELDS:
                        self.logger.critical(
                            (
                                'Improper formatting at classification rejection log file line %i; '
                                'expected %i columns, got %i'
                            ) % (i, REJ_LOG_FIELDS, len(data)),
                        )
                        sys.exit(1)
                    if data[0] != TRANSCRIPT:
                        continue
                    discarded_at_classification.add(data[1])

        self.num_no_chains = len(discarded_at_classification)
        num_tr: str = len(self.all_initial_transcripts) + self.num_no_chains
        self.num_orth = len(self.tr2orths)
        self.fr_orth = to_perc(self.num_orth, num_tr)
        self.num_o2o = sum(len(x) == 1 for x in self.tr2orths.values())
        self.fr_o2o = to_perc(self.num_o2o, num_tr)
        self.num_no_orth = sum(x not in self.tr2orths for x in self.all_initial_transcripts)
        self.fr_no_orth = to_perc(self.num_no_orth, num_tr)
        self.fr_no_chains = to_perc(self.num_no_chains, num_tr)

        ## throw a warning if <50% of transcrips do not have an ortholog
        if self.fr_no_orth > MIN_ORTH_PERCENTAGE:
            warning_message: str = CLASSIFICATION_WARNING.format(
                num=self.num_no_orth, 
                lvl='transcripts',
                fraction=self.fr_no_orth
            )
            log_message = warning_message
            raise_warning = True
        else:
            warning_message: str = ''
            log_message: str = CLASSIFICATION_OK.format(
                num=self.num_no_orth,
                lvl='transcripts',
                fraction=self.fr_no_orth
            )
            raise_warning = False
        elem_trs: ElementaryCheckResult = ElementaryCheckResult(
            is_warning=raise_warning, warning_message=warning_message, log_message=log_message
        )
        messages.append(elem_trs)

        ## parse isoforms file if present
        if self.isoforms is not None:
            self._parse_isoform_file()
            num_gene: int = sum(
                any(y in self.all_initial_transcripts for y in x) for x in self.gene2tr
            )
            self.gene_num_orth: int = sum(
                any(y in self.tr2orths for y in x) for x in self.gene2tr.values()
            )
            self.gene_fr_orth = to_perc(self.gene_num_orth, num_gene)
            self.gene_num_o2o = sum(
                all(self.tr2orths[y] == 1 for y in x) for x in self.gene2tr.values()
            )
            self.gene_fr_o2o = to_perc(self.gene_num_o2o, num_gene)
            ## throw a warning if <50% genes have an ortholog
            if self.gene_fr_no_orth > MIN_ORTH_PERCENTAGE:
                warning_message: str =  CLASSIFICATION_WARNING.format(
                    num=self.gene_num_o2o,
                    lvl='genes',
                    fraction=self.gene_fr_o2o
                )
                log_message: str = warning_message
                raise_warning: bool = True
            else:
                warning_message: str = ''
                log_message += '\n' + CLASSIFICATION_OK.format(
                    num=self.gene_num_o2o,
                    lvl='genes',
                    fraction=self.gene_fr_o2o
                )
                raise_warning: bool = False
            elem_genes: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=raise_warning, warning_message=warning_message, log_message=log_message
            )
            messages.append(elem_genes)

        result: SanityCheckResult = SanityCheckResult(
            step='classification', messages=messages
        )
        return result

    def check_preprocessing(self) -> SanityCheckResult:
        """
        Checks preprocessing rejection log. 
        Throws a warning if >33% of projections having reached this point 
        were discarded before alignment. 
        """
        ## get the list of transcripts that reached this point
        if not os.path.exists(self.memory_report):
            self.logger.critical(
                'Memory report file %s does not exist' % self.memory_report
            )
            sys.exit(1)
        if not os.path.exists(self.prepr_rej_log):
            self.logger.warning(
                (
                    'Preprocessing step rejection log file %s does not exist; '
                    'assuming no items were discarded at the preprocessing step'
                ) % self.prepr_rej_log
            )
            result: SanityCheckResult = SanityCheckResult(
                step='preprocessing',
                messages=[
                    ElementaryCheckResult(
                        is_warning=False, warning_message='', log_message='No items were discarded at preprocessing step'
                    )
                ]
            )
            return result

        accepted_num: int = 0
        with open(self.memory_report, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if data[0] == MEM_FILE_HEADER:
                    continue
                if len(data) != MEM_FILE_FIELDS:
                    self.logger.critical(
                        (
                            'Improper formatting at memory file line %s; '
                            'expected %i fields, got %s'
                        ) % (i, MEM_FILE_FIELDS, len(data))
                    )
                    sys.exit(1)
                accepted_num += 1
        ## get the list of transcripts discarded at this step 
        total_rej_num: int = 0
        rej_num: int = 0
        with open(self.prepr_rej_log, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != REJ_LOG_FIELDS:
                    self.logger.critical(
                        (
                            'Improper formatting at preprocessing step rejection log file line %s; '
                            'expected %i fields, got %s'
                        ) % (i, REJ_LOG_FIELDS, len(data))
                    )
                total_rej_num += 1
                if data[4] == CHAIN_LIMIT_REASON:
                    continue
                rej_num += 1
        total_num: int = accepted_num + rej_num
        rej_perc: float = to_perc(rej_num, total_num)
        if rej_perc > MAX_REJ_AT_PREPR_PERCENTAGE:
            warning_message = PREPROCESSING_WARNING.format(total_num, rej_num, rej_perc)
            log_message = warning_message
            raise_warning = True
        else:
            warning_message: str = ''
            log_message = PREPROCESSING_OK.format(total_num, rej_num, rej_perc)
            raise_warning = False

        result: SanityCheckResult = SanityCheckResult(
            step='preprocessing',
            messages=[
                ElementaryCheckResult(
                    is_warning=raise_warning, warning_message=warning_message, log_message=log_message
                )
            ]
        )

        return result

    def check_loss_summary(self) -> None:
        """
        Checks gene loss summary. Throws a warning if:
        a) <50% of original projections are considered present;
        b) <50% of reference transcripts are considered present;
        c) <50% of reference genes are considered present.
        """
        self._parse_non_orthologous()

        if not os.path.exists(self.loss_summary):
            self.logger.critical('Gene loss summary file %s does not')
            sys.exit(1)
        messages: List[ElementaryCheckResult] = []
        with open(self.loss_summary, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != LOSS_FILE_FIELDS:
                    self.logger.critical(
                        (
                            'Improper formatting at loss summary file line %i; '
                            'expected %i fields, got %i'
                        ) % (i, LOSS_FILE_FIELDS, len(data))
                    )
                    sys.exit(1)
                if data[0] == LOSS_HEADER:
                    continue
                level, name, status = data
                if level == PROJECTION:
                    # self.proj2loss[status] += 1
                    if name in self.non_orthologous:
                        ## ignore Ns, PPgenes, and paralogs for the warning message
                        continue
                    if status == 'N':
                        continue
                    self.proj2loss[status] += 1
                if level == TRANSCRIPT:
                    if status in IGNORED_STATUSES:
                        continue
                    self.tr2loss[status] += 1
                if level == GENE:
                    if status in IGNORED_STATUSES:
                        continue
                    self.gene2loss[name] = status
        ## throw a warning for each loss summary level
        ## 1) projections
        num_proj: int = sum(self.proj2loss.values())
        num_acc_proj: int = sum(y for x,y in self.proj2loss.items() if x in self.intact_classes)
        fr_acc_proj: float = to_perc(num_acc_proj, num_proj)
        if fr_acc_proj < MIN_INTACT_PROJ_PERCENTAGE:
            warning_message: str = LOSS_WARNING.format(
                lvl='projections',
                total=num_proj, 
                num_valid=num_acc_proj, 
                fraction=fr_acc_proj, 
                classes=','.join(self.intact_classes)
            )
            log_message = warning_message
            elem_proj: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True, warning_message=warning_message, log_message=log_message
            )
        else:
            warning_message: str = ''
            log_message = LOSS_OK.format(
                lvl='projections',
                total=num_proj, 
                num_valid=num_acc_proj, 
                fraction=fr_acc_proj, 
                classes=','.join(self.intact_classes)
            )
            elem_proj: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False, warning_message=warning_message, log_message=log_message
            )
        messages.append(elem_proj)
        ## 2) reference transcripts
        num_tr: int = sum(self.tr2loss.values())
        num_acc_tr: int = sum(y for x,y in self.tr2loss.items() if x in self.intact_classes)
        fr_acc_tr: float = to_perc(num_acc_tr, num_tr)
        if fr_acc_tr < MIN_INTACT_PROJ_PERCENTAGE:
            warning_message: str = LOSS_WARNING.format(
                lvl='reference transcripts',
                total=num_tr, 
                num_valid=num_acc_tr, 
                fraction=fr_acc_tr, 
                classes=','.join(self.intact_classes)
            )
            log_message = warning_message
            elem_tr: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True, warning_message=warning_message, log_message=log_message
            )
        else:
            warning_message: str = ''
            log_message = LOSS_OK.format(
                lvl='reference transcripts',
                total=num_tr, 
                num_valid=num_acc_tr, 
                fraction=fr_acc_tr, 
                classes=','.join(self.intact_classes)
            )
            elem_tr: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False, warning_message=warning_message, log_message=log_message
            )
        messages.append(elem_tr)
        ## 3) genes (if isoform file was provided)
        if self.isoforms is not None:
            num_gene: int = sum(self.gene2loss.values())
            num_acc_gene: int = sum(y for x,y in self.gene2loss.items() if x in self.intact_classes)
            fr_acc_gene: int = to_perc(num_acc_gene, num_gene)
            if fr_acc_tr < MIN_INTACT_PROJ_PERCENTAGE:
                warning_message: str = LOSS_WARNING.format(
                    lvl='reference genes',
                    total=num_gene, 
                    num_valid=num_acc_gene, 
                    fraction=fr_acc_gene, 
                    classes=','.join(self.intact_classes)
                )
                log_message = warning_message
                elem_gene: ElementaryCheckResult = ElementaryCheckResult(
                    is_warning=True, warning_message=warning_message, log_message=log_message
                )
            else:
                warning_message: str = ''
                log_message = LOSS_OK.format(
                    lvl='reference genes',
                    total=num_gene, 
                    num_valid=num_acc_gene, 
                    fraction=fr_acc_gene, 
                    classes=','.join(self.intact_classes)
                )
                elem_gene: ElementaryCheckResult = ElementaryCheckResult(
                    is_warning=False, warning_message=warning_message, log_message=log_message
                )
            messages.append(elem_gene)
        
        result: SanityCheckResult = SanityCheckResult(
            step='gene loss summary', messages=messages
        )
        return result

    def check_orthology_resolution(self) -> None:
        """
        Checks orthology resolution results. Throws a warning if:
        a) <50% of genes are classified as one:one;
        b) >=20% of genes are classified as one:zero 
        """
        if not os.path.exists(self.orth_classes):
            self.logger.critical(
                'Orthology summary file %s does not exist' % self.orth_classes
            )
            sys.exit(1)
        with open(self.orth_classes, 'r') as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if data[0] == ORTH_CLASS_HEADER:
                    continue
                if len(data) != ORTH_CLASS_FIELDS:
                    self.logger.critical(
                        (
                            'Improper formatting at orthology summary file line %i; '
                            'expected %i fields, got %i'
                        ) % (i, ORTH_CLASS_FIELDS, len(data))
                    )
                # ref_gene: str = data[0]
                # query_gene: str = data[2]
                status: str = data[4]
                self.orth_status_counter[status] += 1
                ## TODO: For sanity check, do not account for genes discarded at the initial steps
                ## (technical reasons, etc.)
                ## NOTE: those with no orthologous predictions at the initial step should still be accounted for!
        messages: List[ElementaryCheckResult] = []
        num_gene: int = sum(self.orth_status_counter.values())
        num_one2one: int = self.orth_status_counter.get('one2one', 0)
        fr_one2one: float = to_perc(num_one2one, num_gene )
        if fr_one2one < MIN_ONE2ONE_PERCENTAGE:
            warning_message: str = RESOLUTION_WARNING.format(
                total=num_gene,
                num=num_one2one,
                fraction=fr_one2one,
                type='1:1'
            )
            log_message: str = warning_message
            elem_one2one: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True, warning_message=warning_message, log_message=log_message
            )
        else:
            warning_message: str = ''
            log_message: str = RESOLUTION_OK.format(
                total=num_gene,
                num=num_one2one,
                fraction=fr_one2one,
                type='1:1'
            )
            elem_one2one: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False, warning_message=warning_message, log_message=log_message
            )
        messages.append(elem_one2one)
        num_one2zero: int = self.orth_status_counter.get('one2zero', 0)
        fr_one2zero: float = to_perc(num_one2zero, num_gene)
        if fr_one2zero > MAX_ONE2ZERO_PERCENTAGE:
            warning_message: str = RESOLUTION_WARNING.format(
                total=num_gene,
                num=num_one2zero,
                fraction=fr_one2zero,
                type='1:0'
            )
            log_message: str = warning_message
            elem_one2one: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=True, warning_message=warning_message, log_message=log_message
            )
        else:
            warning_message: str = ''
            log_message: str = RESOLUTION_OK.format(
                total=num_gene,
                num=num_one2zero,
                fraction=fr_one2zero,
                type='1:0'
            )
            elem_one2zero: ElementaryCheckResult = ElementaryCheckResult(
                is_warning=False, warning_message=warning_message, log_message=log_message
            )
            messages.append(elem_one2zero)

        result: SanityCheckResult = SanityCheckResult(step='orthology resolution', messages=messages)
        return result

    def summary(
        self, 
        ref_2bit: str, 
        query_2bit: str,
        chains: str,
        ref_annot: str,
        output: str
    ) -> str:
        """Summarizes most crucial TOGA2 result statistics"""
        isoforms_line: str = (
            ISOFORMS_LINE.format(self.isoforms) if self.isoforms is not None 
            else ISOFORMS_CAVEAT
        )
        ## 1) orthology classification
        if self.isoforms is not None:
            gene_level_attrs: List[Union[int, float]] = [
                self.gene_num_orth, self.gene_fr_orth,
                self.gene_num_o2o, self.gene_fr_o2o,
                self.gene_num_no_orth, self.gene_fr_no_orth,
            ]
        else:
            gene_level_attrs: List[Union[int, float]] = [
            self.num_orth, self.fr_orth,
            self.num_o2o, self.fr_o2o,
            self.num_no_orth, self.fr_no_orth,
            ]
        orth_class: str = CLASSIFICATION_BOILERPLATE.format(
            len(self.all_initial_transcripts),
            self.num_orth, self.fr_orth,
            self.num_o2o, self.fr_o2o,
            self.num_no_orth, self.fr_no_orth,
            self.num_no_chains, self.fr_no_chains,
            self.num_ppgene,
            *gene_level_attrs
        )

        ## 2) gene loss summary
        num_proj: str = sum(self.proj2loss.values())
        proj2fr: Dict[str, float] = {x: to_perc(y, num_proj) for x, y in self.proj2loss.items()}
        num_proj_present: int = sum(
            y for x, y in self.proj2loss.items() if x in self.intact_classes
        )
        fr_proj_present: float = to_perc(num_proj_present, num_proj)
        num_proj_missing: int = num_proj - num_proj_present
        fr_proj_missing: float = 100.0 - fr_proj_present
        num_tr: int = sum(self.tr2loss.values())
        tr2fr: Dict[str, float] = {x: to_perc(y, num_tr) for x, y in self.tr2loss.items()}
        num_tr_present: int = sum(
            y for x, y in self.tr2loss.items() if x in self.intact_classes
        )
        fr_tr_present: float = to_perc(num_tr_present, num_tr)
        num_tr_missing: int = num_tr - num_tr_present
        fr_tr_missing: float = 100.0 - fr_tr_present
        if self.isoforms is not None:
            num_gene: str = sum(self.gene2loss.values())
            gene2fr: Dict[str, float] = {x: to_perc(y, num_gene) for x, y in self.gene2loss.items()}
            num_gene_present: int = sum(
                y for x, y in self.gene2loss.items() if x in self.intact_classes
            )
            fr_gene_present: float = to_perc(num_gene_present, num_gene)
            num_gene_missing: int = num_gene - num_gene_present
            fr_gene_missing: float = 100.0 - fr_gene_present
        else:
            gene2fr: Dict[str, float] = tr2fr
            num_gene_present: int = num_tr_present
            fr_gene_present: float = fr_tr_present
            num_gene_missing: int = num_tr_missing
            fr_gene_missing: float = fr_tr_missing
        loss_summary: str = LOSS_SUMMARY_BOILERPLATE.format(
            ','.join(self.intact_classes),
            self.proj2loss.get('FI', 0), proj2fr.get('FI', 0.0),
            self.proj2loss.get('I', 0), proj2fr.get('I', 0.0),
            self.proj2loss.get('PI', 0), proj2fr.get('PI', 0.0),
            self.proj2loss.get('UL', 0), proj2fr.get('UL', 0.0),
            self.proj2loss.get('L', 0), proj2fr.get('L', 0.0),
            self.proj2loss.get('M', 0), proj2fr.get('M', 0.0),
            num_proj_present, fr_proj_present,
            num_proj_missing, fr_proj_missing,
            self.tr2loss.get('FI', 0), tr2fr.get('FI', 0.0),
            self.tr2loss.get('I', 0), tr2fr.get('I', 0.0),
            self.tr2loss.get('PI', 0), tr2fr.get('PI', 0.0),
            self.tr2loss.get('UL', 0), tr2fr.get('UL', 0.0),
            self.tr2loss.get('L', 0), tr2fr.get('L', 0.0),
            self.tr2loss.get('M', 0), tr2fr.get('M', 0.0),
            num_tr_present, fr_tr_present,
            num_tr_missing, fr_tr_missing,
            self.gene2loss.get('FI', 0), gene2fr.get('FI', 0.0),
            self.gene2loss.get('I', 0), gene2fr.get('I', 0.0),
            self.gene2loss.get('PI', 0), gene2fr.get('PI', 0.0),
            self.gene2loss.get('UL', 0), gene2fr.get('UL', 0.0),
            self.gene2loss.get('L', 0), gene2fr.get('L', 0.0),
            self.gene2loss.get('M', 0), gene2fr.get('M', 0.0),
            num_gene_present, fr_gene_present,
            num_gene_missing, fr_gene_missing
        )

        ## 3) orthology resoltion
        num_ref_genes: int = sum(self.orth_status_counter.values())
        query_orth, query_orphan = self._get_query_gene_num()
        num_query_genes: int = query_orth + query_orphan
        fr_orth_query: float = to_perc(query_orth, num_query_genes)
        fr_orphan_query: float = 100.0 - fr_orth_query
        orth2fr: Dict[str, float] = {
            x: to_perc(y, num_ref_genes) for x, y in self.orth_status_counter.items()
        }
        orth_summary: str = ORTHOLOGY_BOILERPLATE.format(
            num_ref_genes, num_query_genes, 
            query_orth, fr_orth_query, query_orphan, fr_orphan_query,
            self.orth_status_counter.get('one2one', 0), orth2fr.get('one2one', 0.0),
            self.orth_status_counter.get('one2many', 0), orth2fr.get('one2many', 0.0),
            self.orth_status_counter.get('many2one', 0), orth2fr.get('many2one', 0.0),
            self.orth_status_counter.get('many2many', 0), orth2fr.get('many2many', 0.0),
            self.orth_status_counter.get('one2zero', 0), orth2fr.get('one2zero', 0.0)
        )

        summary: str = SUMMARY_BOILERPLATE.format(
            br=BREAK_LINE,
            ref_2bit=ref_2bit,
            query_2bit=query_2bit,
            chains=chains,
            ref_annot=ref_annot,
            output=output,
            isoforms_line=isoforms_line,
            orthology_class_report=orth_class,
            loss_summary=loss_summary,
            orthology_res_report=orth_summary
        )

        return summary

