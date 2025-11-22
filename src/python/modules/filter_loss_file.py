#!/usr/bin/env python3

""" """

from typing import List, Optional, Set

import click

from .constants import Headers
from .shared import (
    CONTEXT_SETTINGS,
    CommandLineManager,
    base_proj_name,
    parse_single_column,
)

HEADER: str = "level"
PROJECTION: str = "PROJECTION"


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    "init_loss_summary", type=click.File("r", lazy=True), metavar="INIT_LOSS_SUMMARY"
)
@click.argument(
    "query_bed_file", type=click.File("r", lazy=True), metavar="FINAL_BED_FILE"
)
@click.argument(
    "final_loss_summary", type=click.File("w", lazy=True), metavar="FINAL_LOSS_SUMMARY"
)
@click.option(
    "--paralogs",
    "-p",
    type=click.File("r", lazy=True),
    metavar="PARALOG_FILE",
    default=None,
    show_default=True,
    help="A single-column file containing paralogous projections' identifiers",
)
@click.option(
    "--processed_pseudogenes",
    "-pp",
    type=click.File("r", lazy=True),
    metavar="PROCESSED_PSEUDOGENE_FILE",
    default=None,
    show_default=True,
    help="A single-column file containing processed pseudogene projections' identifiers",
)
@click.option(
    "--log_name",
    "-ln",
    type=str,
    metavar="STR",
    default=None,
    show_default=True,
    help="Logger name to use; relevant only upon main class import",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls execution verbosity",
)
class LossFileFilter(CommandLineManager):
    """
    A straightforward loss summary filtering class. Given a full loss summary file
    and a final annotation BED file, picks loss records only for projections present
    in the annotation. Transcript- and gene-level records are not filtered and reported
    as in the original file.

    \b
    Arguments are:
    \t* INIT_LOSS_SUMMARY is a full loss summary file produced by TOGA2 (`meta/loss_summary_extended.tsv` by default);
    \t* FINAL_BED_FILE is a final query annotation file in BED12 format (`query_annotation.tsv` or `query_annotation.with_utrs.bed`);
    \t* FINAL_LOSS_SUMMARY is a path to the final, filtered loss summary file (`loss_summary.tsv` in TOGA2 by default)
    """

    __slots__ = ("loss_in", "bed_file", "loss_out", "paralogs", "ppgenes")

    def __init__(
        self,
        init_loss_summary: click.File,
        query_bed_file: click.File,
        final_loss_summary: click.File,
        paralogs: Optional[click.File],
        processed_pseudogenes: Optional[click.File],
        log_name: Optional[str],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)
        self.loss_in: click.File = init_loss_summary
        self.bed_file: click.File = query_bed_file
        self.loss_out: click.File = final_loss_summary
        self.paralogs: Set[str] = parse_single_column(paralogs)
        self.ppgenes: Set[str] = parse_single_column(processed_pseudogenes)

        self.run()

    def run(self) -> None:
        """Entry point"""

        projections: Set[str] = set()
        for i, line in enumerate(self.bed_file, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 12:
                self._die(
                    (
                        "Improper final Bed file formatting at line %i; "
                        "expected 12 fields, got %i"
                    )
                    % (i, len(data))
                )
            name: str = base_proj_name(data[3])
            projections.add(name)
        self.loss_out.write(Headers.LOSS_FILE_HEADER)
        for i, line in enumerate(self.loss_in, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 3:
                self._die(
                    (
                        "Improper loss file formatting at line %i; "
                        "expected 3 fields, got %i"
                    )
                    % (i, len(data))
                )
            if data[0] == HEADER:
                continue
            if data[0] != PROJECTION:
                self.loss_out.write(line)
                continue
            if (
                data[1] not in projections
                and base_proj_name(data[1]) not in projections
            ):
                continue
            if data[1] in self.paralogs:
                data[1] += "#paralog"
                line = "\t".join(data) + "\n"
            if data[1] in self.ppgenes:
                data[1] += "#retro"
                line = "\t".join(data) + "\n"
            self.loss_out.write(line)


if __name__ == "__main__":
    LossFileFilter()
