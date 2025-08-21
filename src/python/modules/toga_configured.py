#!/usr/bin/env python3

"""
Wrapper over TogaMain for file-configured runs
"""

from .constants import TOGA2_SLOT2ARG
from .shared import CommandLineManager
from toga_main import __version__, TogaMain
from typing import Optional, Tuple, Union

import click
import os

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__all__ = (None)

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
TOGA2: str = os.path.join(LOCATION, 'toga2.py')
COL_NUM: int = 2
REF: str = 'ref_2bit'
QUERY: str = 'query_2bit'
CHAIN: str = 'chain_file'
ANNOT: str = 'ref_annotation'
MANDATORY_ARGS: Tuple[str] = (REF, QUERY, CHAIN, ANNOT)
TRUE: str = 'True'
FALSE: str = 'False'
NONE: str = 'None'
PROJ_ARG_FILE: str = 'project_args.json'
VERSION: str = 'version'

class Toga2ConfiguredLauncher(CommandLineManager):
    __slots__ = ()

    def __init__(self, config_file: click.File, override: Optional[Union[str, None]]) -> None:
        pass