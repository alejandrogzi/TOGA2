#!/usr/bin/env python3

"""
A TOGA2 test function 
"""

from modules.shared import get_upper_dir
from modules.toga_main import TogaMain

import os

TOGA2_ROOT: str = get_upper_dir(__file__, 4)
TEST_DIR: str = os.path.join(TOGA2_ROOT, 'test_input')
MICRO_REF: str = os.path.join(TEST_DIR, 'hg38.micro_sample.2bit')
MICRO_QUERY: str = os.path.join(TEST_DIR, 'q2bt_micro_sample2.bit')
MICRO_CHAIN: str = os.path.join(TEST_DIR, 'align_micro_sample.chain')
MICRO_ANNOT: str = os.path.join(TEST_DIR, 'annot_micro_sample.bed')


def test_micro(output: str) -> None:
    TogaMain(
        ref_2bit=MICRO_REF,
        query_2bit=MICRO_QUERY,
        chain_file=MICRO_CHAIN,
        ref_annotation=MICRO_ANNOT,
        output=output
    )
