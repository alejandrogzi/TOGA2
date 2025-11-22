#!/usr/bin/env python3

"""
Integrates TOGA2 results obtained with multiple references for a single query
to produce a combined, multi-reference TOGA2 annotation
"""

import gzip
import json
import os
from collections import defaultdict
from dataclasses import dataclass
from shutil import which
from typing import Dict, Iterable, List, Optional, Set, TextIO, Union

import networkx as nx

from .cesar_wrapper_constants import CLASS_TO_COL, CLASS_TO_NUM
from .constants import Constants
from .shared import (
    CommandLineManager,
    base_proj_name,
    get_proj2trans,
    get_upper_dir,
    intersection,
)

TOGA2_ROOT: str = get_upper_dir(__file__, 4)
BIN: str = os.path.join(TOGA2_ROOT, "bin")
DEFAULT_BIGBED2BED: str = os.path.join(BIN, "bigBedToBed")
LOCATION: str = os.path.dirname(os.path.abspath(__file__))
MAKE_IX_SCRIPT: str = os.path.join(LOCATION, "get_names_from_bed.py")

ALL: str = "ALL"
BED_FIELD_NUM: int = 12
COL2CLASS: Dict[str, str] = {y: x for x, y in CLASS_TO_COL.items()}
EXON_HEADER: str = "projection"
INTACT_EXON: str = "I"
SUPPORTED: str = "CHAIN_SUPPORTED"
PLUS: str = "+"

## check the NetworkX version
nx_v: str = nx.__version__
v_split: List[str] = [x for x in nx_v.split(".") if x.isnumeric()]
if len(v_split) > 1:
    NX_VERSION: float = float(f"{v_split[0]}.{v_split[1]}")
else:
    NX_VERSION: float = float(v_split[0])


class ReferenceBundle:  ## initialized from a JSON object
    """
    Data class for reference- and query annotation data
    bound to a certain TOGA2 reference
    """

    __slots__ = (
        "query_bed",
        "exon_meta",
        "ref_isoforms",
        "paralogs",
        "ppgenes",
        "ucsc_bigbed",
        "protein_file",
        "nucleotide_file",
        "priority",
    )

    def __init__(
        self,
        kwargs: Dict[str, Union[str, os.PathLike, int]],
        priority: Optional[int] = None,
    ) -> None:
        self.query_bed: str = kwargs["query_bed"]
        self.exon_meta: str = kwargs["exon_meta"]
        self.ref_isoforms: Union[str, None] = kwargs.get("reference_isoforms", None)
        self.paralogs: Union[str, None] = kwargs.get("paralog_list", None)
        self.ppgenes: Union[str, None] = kwargs.get("processed_pseudogene_list", None)
        self.ucsc_bigbed: Union[str, None] = kwargs.get("ucsc_bigbed", None)
        self.protein_file: Union[str, None] = kwargs.get("protein_file", None)
        self.nucleotide_file: Union[str, None] = kwargs.get("nucleotide_file", None)
        self.priority: int = kwargs.get("priority", priority)


@dataclass
class BedRecord:
    """
    Data class containing data for a single BED file record
    """

    __slots__ = (
        "name",
        "ref",
        "chrom",
        "start",
        "end",
        "strand",
        "loss_status",
        "lines",
        "exons",
    )

    def __init__(
        self,
        name: str,
        ref: str,
        chrom: str,
        start: int,
        end: int,
        strand: bool,
        loss_status: str,
        lines: List[str],
    ) -> None:
        self.name: str = name
        self.ref: str = ref
        self.chrom: str = chrom
        self.start: int = start
        self.end: int = end
        self.strand: bool = strand
        self.loss_status: str = loss_status
        self.lines: List[str] = lines
        self.exons: List[str] = []

    def return_bed_line(self, prefix: Union[str] = "") -> Iterable[str]:
        """Returns the initial BED line for the projection"""
        for line in self.lines:
            yield line.format(f"{prefix}.{self.name}", CLASS_TO_COL[self.loss_status])


@dataclass
class ExonRecord:
    """Exon coordinate record"""

    ## TODO: Salvage a similar class from preprocessing code
    __slots__ = ("projection", "num", "chrom", "start", "end", "strand")

    projection: str
    num: int
    chrom: str
    start: int
    end: int
    strand: bool


class AnnotationIntegrator(CommandLineManager):
    __slots__ = (
        "ref_data",
        "query_projections",
        "query_proj2ref",
        "query_annotation",
        "ref_proj2gene",
        "paralog_pool",
        "ppgene_pool",
        "graph",
        "discarded_items",
        "final_projections",
        "accepted_statuses",
        "output",
        "gene_tsv",
        "gene_bed",
        "projection_bed",
        "protein_file",
        "nucleotide_file",
        "ucsc_dir",
        "bigbed_stub",
        "bigbed",
        "ix",
        "ixx",
        "prefix",
        "skip_ucsc",
        "has_ucsc_data",
        "bigbedtobed_binary",
        "bedtobigbed_binary",
        "ixixx_binary",
        "schema",
        "chrom_sizes",
        "bed_index",
    )

    def __init__(
        self,
        ref_data: Union[str, os.PathLike],
        output: Union[str, os.PathLike],
        accepted_statuses: str,
        prefix: str,
        skip_ucsc: bool,
        chrom_sizes: Union[str, os.PathLike, None],
        bigbedtobed_binary: Union[str, os.PathLike],
        bedtobigbed_binary: Union[str, os.PathLike],
        ixixx_binary: Union[str, os.PathLike],
        verbose: Optional[bool] = False,
    ) -> None:
        self.v: bool = verbose
        self.set_logging()

        ## parse the input dictionary
        self.ref_data: Dict[str, ReferenceBundle] = {}
        with open(ref_data, "r") as h:
            bundles: Dict[str, Dict[str, str]] = json.loads(h.read())
            for i, (species, data) in enumerate(bundles.items()):
                self.ref_data[species] = ReferenceBundle(data, priority=i)
        if accepted_statuses == ALL:
            self.accepted_statuses: List[str] = Constants.ALL_LOSS_SYMBOLS
        else:
            self.accepted_statuses: List[str] = [
                x for x in accepted_statuses.split(",") if x
            ]  ## TODO: Add sanity checks
        self.query_projections: Dict[str, BedRecord] = {}
        self.query_proj2ref: Dict[str, str] = {}
        self.query_annotation: Dict[str, List[str]] = defaultdict(list)
        self.ref_proj2gene: Dict[str, str] = {}
        self.paralog_pool: Set[str] = set()
        self.ppgene_pool: Set[str] = set()
        self.graph: nx.Graph = nx.Graph()
        self.discarded_items: Set[str] = set()
        self.final_projections: Set[str] = set()

        self.skip_ucsc: bool = skip_ucsc
        self.bigbedtobed_binary: Union[str, os.PathLike, None] = bigbedtobed_binary
        self.bedtobigbed_binary: Union[str, os.PathLike, None] = bedtobigbed_binary
        self.ixixx_binary: Union[str, os.PathLike, None] = ixixx_binary
        self.chrom_sizes: Union[str, os.PathLike, None] = chrom_sizes
        if not skip_ucsc:
            self._check_binaries()
        if self.chrom_sizes is None and not (
            self.skip_ucsc or all(x.ucsc_bigbed is None for x in self.ref_data.values())
        ):
            self._die("Query chromosome size file was not provided")
        self.schema: str = os.path.join(TOGA2_ROOT, "supply", "bb_schema_toga2.as")
        self.has_ucsc_data: List[str] = []

        self.output: str = output
        self.gene_tsv: str = os.path.join(output, "query_genes.tsv")
        self.gene_bed: str = os.path.join(output, "query_genes.bed")
        self.projection_bed: str = os.path.join(output, "query_annotation.bed")
        self.ucsc_dir: str = os.path.join(output, "ucsc_browser_files")

        self.protein_file: str = os.path.join(output, "protein.fa")
        self.nucleotide_file: str = os.path.join(output, "nucleotide.fa")

        self.bigbed_stub: str = os.path.join(self.ucsc_dir, f"{prefix}.bed")
        self.bigbed: str = os.path.join(self.ucsc_dir, f"{prefix}.bb")
        self.ix: str = os.path.join(self.ucsc_dir, f"{prefix}.ix")
        self.ixx: str = os.path.join(self.ucsc_dir, f"{prefix}.ixx")
        self.bed_index: str = os.path.join(self.ucsc_dir, f"{prefix}.ix.txt")

    def run(self) -> None:
        """Main execution method"""
        self._mkdir(self.output)
        self.process()
        self.infer_genes()
        self.pick_best_isoforms()
        if not self.skip_ucsc:
            self._mkdir(self.ucsc_dir)
        self.prepare_ucsc_file()
        self.prepare_sequence_files()

    def process(self) -> None:
        """Extract data for every reference directory"""
        for species in self.ref_data:
            self.read_annotation(species)
            self.read_exon_meta(species)
            self.read_ref_isoforms(species)
            self.read_paralogs(species)
            self.read_ppgenes(species)

    def _check_binaries(self) -> None:
        """
        Checks the provided UCSC binaries,
        looks for an alternative if needed
        """
        for attr, default_name in Constants.BINARIES_TO_CHECK.items():
            if attr not in self.__slots__:
                continue
            binary: Union[str, None] = getattr(self, attr)
            if binary is not None:
                self._to_log("Testing %s binary at %s" % (default_name, binary))
                if os.access(binary, os.X_OK):
                    self._to_log(
                        "The provided binary is executable; using the stated %s instance"
                        % default_name
                    )
                    continue
                else:
                    self._to_log(
                        (
                            "%s binary at %s does not seem to be executable; "
                            "looking for alternatives"
                        )
                        % (default_name, binary),
                        "warning",
                    )
            else:
                self._to_log(
                    (
                        "No %s executable was not provided; "
                        "looking for alternatives in $PATH"
                    )
                    % default_name
                )
            self._to_log(
                "Checking the default %s executable at TOGA2/bin" % default_name
            )
            default_path: str = os.path.join(BIN, default_name)
            if os.path.exists(default_path):
                self._to_log(
                    "Found default %s executable at %s" % (default_name, default_path)
                )
                if os.access(default_path, os.X_OK):
                    self._to_log(
                        "The default UCSC %s binary is executable; using the default version"
                        % default_name
                    )
                    self.__setattr__(attr, default_path)
                    continue
                else:
                    self._to_log(
                        (
                            "Default UCSC %s binary at %s does not seem to be executable; "
                            "looking for alternatives in $PATH"
                        )
                        % (default_name, default_path),
                        "warning",
                    )
            binary_in_path: Union[str, None] = which(default_name)
            if binary_in_path is not None:
                self._to_log(
                    (
                        "Found bigBedTobed instance at %s; "
                        "checking the execution permissions"
                    )
                    % binary_in_path
                )
                if os.access(binary_in_path, os.X_OK):
                    self._to_log(
                        "The found binary is executable; using the found %s instance"
                        % default_name
                    )
                    self.__setattr__(attr, binary_in_path)
                    continue
                self._die(
                    (
                        "The %s binary found in $PATH at %s is not executable; "
                        "check your $PATH or provide a valid bigBedToBed instance"
                    )
                    % (default_name, binary_in_path)
                )
            self._die(
                (
                    "No %s binary found in $PATH; "
                    "check your $PATH or provide a valid bigBedToBed instance"
                )
                % default_name
            )

    def read_annotation(self, species: str) -> None:
        """
        Parses the reference BED file, filtering the records by their loss status
        """
        self._to_log("Reading annotation file for reference %s" % species)
        file: str = self.ref_data[species].query_bed
        if not os.path.exists(file):
            self._die("File %s does not exist" % file)
        with open(file, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) != BED_FIELD_NUM:
                    self._die(
                        (
                            "Improper formatting at BED file %s line %i; "
                            "expected %i fields, got %i"
                        )
                        % (file, i, BED_FIELD_NUM, len(data))
                    )
                color: str = data[8]
                status: str = COL2CLASS.get(color)
                if status is None:
                    self._die(
                        ("Invalid loss status color at BED file %s line %i: %s")
                        % (file, i, color)
                    )
                # if status not in self.accepted_statuses:
                #     continue
                name: str = data[3]
                chrom: str = data[0]
                start: int = int(data[6])
                end: int = int(data[7])
                strand: bool = data[5] == "+"
                line_template: str = "\t".join(
                    [*data[:3], "{}", *data[4:8], "{}", *data[9:]]
                )
                if name in self.query_projections:
                    if "," in name:
                        self.query_projections[name].lines.append(line_template)
                    else:
                        self._die(
                            "Duplicated non-fragmented entry for reference %s: %s"
                        ) % (species, name)
                else:
                    record: BedRecord = BedRecord(
                        name,
                        species,
                        chrom,
                        start,
                        end,
                        strand,
                        status,
                        [line_template],
                    )
                    self.query_projections[name] = record
                self.query_proj2ref[name] = species
                self.query_annotation[chrom].append(name)

    def read_exon_meta(self, species) -> None:
        """
        Extracts exon coordinates from the exon metadata file for a given species,
        removing missing/deleted and extrapolated exons
        """
        self._to_log("Reading exon metadata file for reference %s" % species)
        file: str = self.ref_data[species].exon_meta
        gzipped: bool = file.endswith(".gz")
        if not os.path.exists(file):
            self._die("File %s does not exist" % file)
        with gzip.open(file, "rb") if gzipped else open(file, "r") as h:
            for line in h:
                if gzipped:
                    line: str = line.decode("utf8")
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                ## TODO: Array length check
                if data[0] == EXON_HEADER:
                    continue
                proj: str = data[0]
                ## do not include projections from transcripts outside of the final annotation files
                if proj not in self.query_projections:
                    continue
                status: str = data[7]
                ## do not account for missing and deleted exons
                if status != INTACT_EXON:
                    continue
                support: str = data[19]
                ## do not account for extrapolated exons
                if support != SUPPORTED:
                    continue
                chrom: str = data[3]
                start: int = int(data[4])
                end: int = int(data[5])
                strand: bool = data[6] == PLUS
                exon_num: int = int(data[1])
                record: str = ExonRecord(proj, exon_num, chrom, start, end, strand)
                self.query_projections[proj].exons.append(record)

    def read_ref_isoforms(self, species: str) -> None:
        """Reads the transcript-to-gene mapping for a given reference"""
        file: str = self.ref_data[species].ref_isoforms
        if file is None:
            self._to_log(
                "No reference isoform file provided for reference %s" % species
            )
            return
        if not os.path.exists(file):
            self._to_log("File %s does not exist; skipping" % file)
            return
        with open(file, "r") as h:
            for line in h:
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                gene, tr = data
                self.ref_proj2gene[tr] = gene

    def read_paralogs(self, species: str) -> None:
        """Extracts paralogous projections' names"""
        file: Union[str, None] = self.ref_data[species].paralogs
        if file is None:
            self._to_log(
                "No paralogous projection list provided for reference %s" % species
            )
            return
        if not os.path.exists(file):
            self._to_log("File %s does not exist; skipping" % file, "warning")
            return
        self._to_log("Reading paralogous projection list for reference %s" % species)
        with open(file, "r") as h:
            for line in h:
                line = line.strip()
                self.paralog_pool.add(line)

    def read_ppgenes(self, species: str) -> None:
        """Extracts processed pseudogene/retrogene projections' names"""
        file: Union[str, None] = self.ref_data[species].ppgenes
        if file is None:
            self._to_log(
                "No retrogene projection list provided for reference %s" % species
            )
            return
        if not os.path.exists(file):
            self._to_log("File %s does not exist; skipping" % file, "warning")
            return
        self._to_log("Reading retrogene projection list for reference %s" % species)
        with open(file, "r") as h:
            for line in h:
                line = line.strip()
                self.ppgene_pool.add(line)

    def infer_genes(self) -> None:
        """
        Infers query genes from combined reference annotations.
        The logic is effectively borrowed from infer_query_genes.py :
        any two projections overlapping by at least one coding base
        """
        self._to_log(
            "Inferring query genes from the combined annotation; might take a while"
        )
        ## create the graph structure
        ## iterate over chromosomes
        for chrom, projections in self.query_annotation.items():
            ## sort the projections
            projections.sort(
                key=lambda x: (
                    self.query_projections[x].start,
                    self.query_projections[x].end,
                )
            )
            ## start iterating over projections
            for i, name_out in enumerate(projections):
                proj_out: BedRecord = self.query_projections[name_out]
                out_start: int = proj_out.start
                out_end: int = proj_out.end
                out_paralog: bool = name_out in self.paralog_pool
                out_ppgene: bool = name_out in self.ppgene_pool
                out_ortholog: bool = not (out_paralog or out_ppgene)
                discarded: bool = False
                edges: List[str] = []
                for j, name_in in enumerate(projections[i + 1 :], start=i):
                    proj_in: BedRecord = self.query_projections[name_in]
                    ## if an item has been already discarded, skip it
                    if name_in in self.discarded_items:
                        continue
                    in_start: int = proj_in.start
                    in_end: int = proj_in.end
                    ## as long as transcripts are sorted properly,
                    ## it's safe to break the inner loop here
                    if in_start >= out_end:
                        break
                    ## if proj_in lies upstream, it must have been already processed
                    if in_end <= out_start:
                        continue
                    ## ignore projections encoded by different strands
                    if proj_in.strand != proj_out.strand:
                        continue
                    ## ignore projections encoded by different strands
                    in_paralog: bool = name_in in self.paralog_pool
                    in_ppgene: bool = name_in in self.ppgene_pool
                    in_ortholog: bool = not (in_paralog or in_ppgene)
                    has_intersection: bool = False
                    for ex_out in proj_out.exons:
                        for ex_in in proj_in.exons:
                            if (
                                intersection(
                                    ex_out.start, ex_out.end, ex_in.start, ex_in.end
                                )
                                > 0
                            ):
                                has_intersection = True
                                break
                        if has_intersection:
                            break
                    if has_intersection:
                        ## ortholog + paralog/pp: discard the non-orthologous prediction
                        if out_ortholog and not in_ortholog:
                            self.discarded_items.add(name_in)
                            continue
                        if not out_ortholog and in_ortholog:
                            self.discarded_items.add(name_out)
                            discarded = True
                            break
                        ## ppgene + paralog: discard the paralog - TODO: Discuss
                        if out_ppgene and in_paralog:
                            self.discarded_items.add(name_in)
                            continue
                        if out_paralog and in_ppgene:
                            self.discarded_items.add(name_out)
                            discarded = True
                            break
                        ## namesakes: add the edge
                        edges.append(name_in)
                if discarded:
                    continue
                ## add the edges to the graph
                if name_out not in self.graph.nodes:
                    self.graph.add_node(name_out)
                for name_in in edges:
                    self.graph.add_edge(name_out, name_in)

    def pick_best_isoforms(self) -> None:
        """
        For each query gene (=overlap graph clique), pick the representatives:
        *   If at least one isoform corresponds to the user-defined loss statuses,
            ignore all other isoforms;
            *   If two or more isoforms are completely identical in terms of coordinates,
                pick the one with the better loss status;
            *   If two or more isoforms have identical layout and loss status, pick the one
                corresponding to a higher priority reference/run;
            *   Otherwise, pick the first encountered representative for this layout-status-priority combination
        *   If no projections with the requested status have been encountered, pick the one
        """
        if NX_VERSION < 2.4:
            components = list(nx.connected_component_subgraphs(self.graph))
        else:
            components = [
                self.graph.subgraph(c) for c in nx.connected_components(self.graph)
            ]
        with (
            open(self.gene_tsv, "w") as gt,
            open(self.gene_bed, "w") as gb,
            open(self.projection_bed, "w") as qb,
        ):
            for component in components:
                ## initialize a temporary storage for initial candidates
                selected: Dict[str, str] = {}
                ## define the best class
                best_status: int = 0
                ## set a semaphore for whether the user-defined loss classes
                ## have been encountered in the clique
                allowed_class_found: bool = False
                for name in component:
                    proj: BedRecord = self.query_projections[name]
                    status: int = CLASS_TO_NUM[proj.loss_status]
                    allowed_status: bool = proj.loss_status in self.accepted_statuses
                    if allowed_status:  ## let this projection in
                        if not allowed_class_found:
                            ## this is the first representative of the allowed loss classes;
                            ## wipe out the previous instances
                            selected.clear()
                        allowed_class_found = True
                    else:
                        ## do not let the items worse than the best currentl
                        if status < best_status:
                            continue
                    ## at this point, this is a likely candidate
                    ## however, chances are an item in exactly the same coordinates
                    ## has been already found
                    if proj.lines[0] in selected:
                        prev_name: str = selected[proj.lines[0]]
                        prev_proj: BedRecord = self.query_projections[prev_name]
                        prev_status: int = CLASS_TO_NUM[prev_proj.loss_status]
                        ## pick the one with the best loss status
                        if prev_status > status:
                            continue
                        ## if it is a tie, go for the more preferred reference
                        elif prev_status == status:
                            species: str = self.query_proj2ref[name]
                            priority: int = self.ref_data[species].priority
                            prev_species: str = self.query_proj2ref[prev_name]
                            prev_priority: int = self.ref_data[prev_species].priority
                            ## if the previous prediction's priority is higher (lower) or equal, keep it
                            if prev_priority <= priority:
                                continue
                        ## otherwise, the new prediction is the winner
                    best_status = max(best_status, status)
                    selected[proj.lines[0]] = name
                ## all the projections have been processed; name the gene and define its coordinates
                all_projs: List[BedRecord] = [
                    self.query_projections[x] for x in selected.values()
                ]
                ## define the coordinates; that's the easy part
                chrom: str = all_projs[0].chrom
                strand: str = "+" if all_projs[0].strand else "-"
                start: int = min([x.start for x in all_projs])
                end: int = max([x.end for x in all_projs])
                ## define the name; a slightly trickier part
                ## first, define how many gene participate in the locus annotation
                genes: Set[str] = {
                    self.ref_proj2gene.get(get_proj2trans(x.name)[0], x.name)
                    for x in all_projs
                }
                ## second, define prefix
                if any(base_proj_name(x.name) in self.paralog_pool for x in all_projs):
                    prefix: str = "paralog_"
                elif any(base_proj_name(x.name) in self.ppgene_pool for x in all_projs):
                    prefix: str = "retro_"
                elif not allowed_class_found:
                    if best_status > CLASS_TO_NUM["M"]:
                        prefix: str = "lost_"
                    else:
                        prefix: str = "missing_"
                else:
                    prefix: str = ""
                ## define the main name; for simplicity, do not bother with chain numbers
                if len(genes) == 1:  ## single gene; assign its name to query locus
                    gene: str = genes.pop()
                elif (
                    len(genes) < 4
                ):  ## up to three genes; combine their names separated by comma
                    gene: str = ",".join(genes)
                else:  ## more than three genes; gene+
                    gene: str = genes.pop() + "+"
                name: str = prefix + gene
                ## done! now, write the output
                ## for gene BED, it's a single line
                gene_bed: str = "\t".join(
                    map(str, [chrom, start, end, name, 0, strand])
                )
                gb.write(gene_bed + "\n")
                ## then, iterate over projection
                for proj in all_projs:
                    basename: str = base_proj_name(proj.name)
                    self.final_projections.add(basename)
                    ## get the reference name and add it as a prefix
                    species: str = self.query_proj2ref[proj.name]
                    ## gene isoform file; write the gene and projections names
                    gt.write(name + "\t" + f"{species}.{proj.name}" + "\n")
                    ## projection BED file; write the original BED line
                    for proj_bed in proj.return_bed_line(prefix=species):
                        qb.write(proj_bed + "\n")

    def prepare_ucsc_file(self) -> None:
        """
        Parses the UCSC BigBed file, extracting the lines for integrated annotation BigBed
        as well as sequences for final sequence files
        """
        if self.skip_ucsc:
            self._to_log("Skipping the joint UCSC BigBed file creation as suggested")
            return
        nuc_seqs: List[str] = []
        prot_seqs: List[str] = []
        longest_word: int = 0
        with open(self.bigbed_stub, "w") as h:
            for ref, data in self.ref_data.items():
                seqs_found: Set[str] = set()
                file: str = data.ucsc_bigbed
                if file is None:
                    self._to_log(
                        "No UCSC BigBed file provided for reference %s; skipping" % ref,
                        "warning",
                    )
                    continue
                if not os.path.exists(file):
                    self._die("UCSC BigBed file %s does not exist" % file)
                self.has_ucsc_data.append(ref)
                read_cmd: str = f"{self.bigbedtobed_binary} {file} stdout"
                output: TextIO = self._exec(
                    read_cmd, "Attempt to read BigBed file %s failed: " % file
                )
                for line in output.split("\n"):
                    # line = line.decode('utf8')
                    data: List[str] = line.strip().split("\t")
                    if not data or not data[0]:
                        continue
                    name: str = data[3]
                    basename: str = base_proj_name(name)
                    if basename not in self.final_projections:
                        continue
                    h.write(line + "\n")
                    longest_word = max(longest_word, len(name))
                    if basename in seqs_found:
                        continue
                    nuc_seq: str = data[34]
                    nuc_seqs.append(f">{ref}.{name}\n{nuc_seq}")
                    prot_seq: str = data[35]
                    prot_seqs.append(f">{ref}.{name}\n{prot_seq}")
                    seqs_found.add(basename)
        ## sort the resulting file, convert it into BigBed format
        bb_cmd: str = (
            f"sort -k1,1 -k2,2n -o {self.bigbed_stub} {self.bigbed_stub} && "
            f"{self.bedtobigbed_binary} -type=bed12+26 {self.bigbed_stub} "
            f"{self.chrom_sizes} {self.bigbed} "
            f"-tab -extraIndex=name -as={self.schema}"
        )
        _ = self._exec(bb_cmd, "BigBed generation failed: ")
        self._to_log("BigBed file successfully created")
        bed_ix_cmd: str = (
            f"{MAKE_IX_SCRIPT} {self.bigbed_stub} | sort -u > {self.bed_index}"
        )
        _ = self._exec(bed_ix_cmd, "BED file indexing failed")
        bigbed_ix_cmd: str = (
            f"{self.ixixx_binary} {self.bed_index} {self.ix} {self.ixx} "
            f"-maxWordLength={longest_word}"
        )
        _ = self._exec(bigbed_ix_cmd, "Index creation for UCSC BigBed file failed: ")
        self._to_log("BigBed index successfully created; removing the temporary files")
        self._rm(self.bigbed_stub)
        self._rm(self.bed_index)
        if prot_seqs:
            with open(self.protein_file, "w") as h:
                for prot in prot_seqs:
                    h.write(prot + "\n")
        if nuc_seqs:
            with open(self.nucleotide_file, "w") as h:
                for nuc in nuc_seqs:
                    h.write(nuc + "\n")

    def prepare_sequence_files(self) -> None:
        """p"""
        seq_files_recorded: bool = bool(self.has_ucsc_data)
        for ref, data in self.ref_data.items():
            if ref in self.has_ucsc_data:
                continue
            prot_file: str = data.protein_file
            if prot_file is None or not os.path.exists(prot_file):
                self._die(
                    (
                        "Reference %s has neither valid UCSC BigBed file "
                        "nor valid protein sequence file"
                    )
                    % ref
                )
            nuc_file: str = data.nucleotide_file
            if nuc_file is None or not os.path.exists(nuc_file):
                self._die(
                    (
                        "Reference %s has neither valid UCSC BigBed file "
                        "nor valid nucleotide sequence file"
                    )
                    % ref
                )
            prot_gzipped: bool = prot_file.endswith(".gz")
            with (
                gzip.open(prot_file, "rb") if prot_gzipped else open(prot_file, "r")
            ) as h:
                prot_seqs: List[str] = self._read_fasta(h, ref)
            nuc_gzipped: bool = nuc_file.endswith(".gz")
            with gzip.open(nuc_file, "rb") if nuc_gzipped else open(nuc_file, "r") as h:
                nuc_seqs: List[str] = self._read_fasta(h, ref)
            write_mode: str = "a" if seq_files_recorded else "w"
            with open(self.protein_file, write_mode) as h:
                for entry in prot_seqs:
                    h.write(entry + "\n")
            with open(self.nucleotide_file, write_mode) as h:
                for entry in nuc_seqs:
                    h.write(entry + "\n")
            seq_files_recorded = True

    def _read_fasta(self, handle: TextIO, ref: str) -> List[str]:
        """
        Reads the sequence file, return the list of FASTA entries corresponding to integrated annotation
        """
        header: str = ""
        seq: str = ""
        sequences: List[str] = []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entry: str = f">{ref}.{header}\n{seq}"
                    sequences.append(entry)
                header = ""
                seq = ""
                proj: str = base_proj_name(header.split()[0])
                if proj not in self.final_projections:
                    continue
                header = proj
            if header:
                seq += line
        if seq:
            entry: str = f">{ref}.{header}\n{seq}"
            sequences.append(entry)

    def set_logging(self) -> None:
        """Sets logging and disables logging propagation"""
        super().set_logging()
        self.logger.propagate = False
