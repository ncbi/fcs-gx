#!/usr/bin/env python3
"""
-----------------------------------------------------------------------------
                             PUBLIC DOMAIN NOTICE
                 National Center for Biotechnology Information

  This software is a "United States Government Work" under the terms of the
  United States Copyright Act.  It was written as part of the author's official
  duties as a United States Government employees and thus cannot be copyrighted.
  This software is freely available to the public for use. The National Library
  of Medicine and the U.S. Government have not placed any restriction on its use
  or reproduction.

  Although all reasonable efforts have been taken to ensure the accuracy and
  reliability of this software, the NLM and the U.S. Government do not and
  cannot warrant the performance or results that may be obtained by using this
  software. The NLM and the U.S. Government disclaim all warranties, expressed
  or implied, including warranties of performance, merchantability or fitness
  for any particular purpose.

  Please cite NCBI in any work or product based on this material.

-----------------------------------------------------------------------------
"""
# pylint: disable=C0301,C0302,C0114,C0103,C0116,C0115,R0914,R0916,R0911,R0902,R1702,R0912,R0913,R0915,W0603
# fmt: off


from __future__ import print_function
from enum import Enum

import sys
import os
import json
import argparse
import re
from collections import namedtuple, defaultdict
from dataclasses import dataclass  # requires python 3.7, or 3.6 with dataclasses
from typing import List, Dict, Optional, Tuple
from math import sqrt
import itertools

TAXONOMY_ANALYSIS_VER_MIN = 1
TAXONOMY_ANALYSIS_VER_MAX = 3
FCS_GENOME_REPORT_VER = 2

global_debug = True
global_debug_seq_id = os.environ.get("GX_ACTION_REPORT_DEBUG_SEQ_ID")
global_discard_small_primaries = not os.environ.get("GX_ACTION_REPORT_DONT_DISCARD_SMALL_PRIMARIES")

# AA: TODO: avoid defining variables in terms of negatives,
# E.g. instead of GX_ACTION_REPORT_NO_FIND_INTEGRANTS=0
#                 GX_ACTION_REPORT_FIND_INTEGRANTS=1

# If this env variable is set, the top organism of the span is calculated from all records;
# if not (default) - only from contam records.
global_top_contam_organism = not os.environ.get("GX_ACTION_REPORT_NO_TOP_CONTAM_ORG")
global_find_integrants     = not os.environ.get("GX_ACTION_REPORT_NO_FIND_INTEGRANTS")

# AA: I think this is in percent.
global_low_coverage_threshold = int(os.environ.get("GX_ACTION_REPORT_LOW_COVERAGE_THRESHOLD", "10"))

# AA: Based on other code, I surmise "PA" in this context means "prok or arch"; I think this is percent.
global_pa_same_kingdom_threshold = int(os.environ.get("GX_ACTION_REPORT_PA_SAME_KINGDOM_THRESHOLD", "1"))

CONTAM_LOW_COV = "contaminant(low-coverage)"

# For tests, as global variables can't be set from importing modules

def set_debug(debug):
    global global_debug
    global_debug = debug


def set_top_contam_organism(val):
    global global_top_contam_organism
    global_top_contam_organism = val


# print to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# debug print to stderr
def dprint(*args, **kwargs):
    if global_debug:
        print(*args, file=sys.stderr, **kwargs)


# print to stderr for global debug or seq_id debug
def dsprint(seq_id, *args, **kwargs):
    if global_debug or (
        global_debug_seq_id
        and (
            seq_id == global_debug_seq_id
            or seq_id.startswith(global_debug_seq_id+"~")
        )
    ):
        print(*args, file=sys.stderr, **kwargs)


def is_same_kingdom(da: str, db: str) -> bool:
    return da[:5] == db[:5]


def is_euk(div) -> bool:
    return any(div.startswith(kdm) for kdm in ("anml", "fung", "plnt", "prst"))


def is_prok(div) -> bool:
    return div.startswith("prok")


def is_prok_or_arch(div) -> bool:
    return div.startswith("prok") or div.startswith("arch")


def is_virus(s) -> bool:
    return "virus" in s


def is_prok_virus(s) -> bool:
    # GP-33245
    # After div name transition second variant should be removed
    return s in ("virs:prokaryotic viruses", "prok|virus")


def is_same_kingdom_pa(da: str, db: str) -> bool:
    return (is_prok_or_arch(da) and is_prok_or_arch(db)) or is_same_kingdom(da, db)


def can_host_virus(div, vir_div) -> bool:
    return div[:5] != "virs:" and vir_div.startswith("virs:prokaryotic" if is_prok(div) else "virs:eukaryotic")


def is_mostly_review_contam(records: List) -> bool:

    def get_len(cond):
        return sum(r.end_pos - r.start_pos + 1 for r in records if is_contam(r.result) and cond(r))

    return (
        get_len(lambda r: r.result == CONTAM_LOW_COV) >=
        get_len(lambda r: r.result != CONTAM_LOW_COV)
    )


#############################################################################
@dataclass
class Taxon:
    tax_id      : int = 0  # e.g. 30033
    div         : str = ""  # tax-div, e.g. flies
    cvg_by_div  : int = 0  # absolute coverage (length) by all alignments in the same div
    cvg_by_tax  : int = 0  # absolute coverage (length) by alignments of this taxon
    score       : int = 0  # aggregate score of the alignments of this taxon.

    # geometric mean of cvg_by_div and cvg_by_tax
    def gmean_cvg(self):
        return sqrt(self.cvg_by_div * self.cvg_by_tax)

    def is_prok_div(self) -> bool:
        return is_prok(self.div)

    # i is starting position for a group of tax-specific columns
    @staticmethod
    def from_row(row: List[str], i: int):
        assert i in [6, 12, 18, 24]
        assert row[i + 1] != ""

        return Taxon(
            int(row[i + 0]),
            row[i + 1],
            int(row[i + 2]),
            int(row[i + 3]),
            int(row[i + 4]),
        )


#############################################################################
# Groups defined in GP-34548, should track changes in classify_taxonomy
result_labels_primary = { "primary-div", "primary-div(virus)", "transposon", }

result_labels_contam = {
    "contaminant(human)", 
    "contaminant(div)",
    "same-kingdom-chimeric",  # AA: apparently intentional; see Slack re: GCA_918807975.1 same-kingdom-chimeric
    "contaminant(virus)",
    "contaminant(cross-kingdom)",
    "contaminant(cross-div)",
}

result_labels_inconclusive = {
    "repeat", 
    "low-coverage",
    "inconclusive",
    "inconclusive(metagenome)",
}

result_labels_valid = result_labels_primary | result_labels_contam | result_labels_inconclusive

# Add local label for contaminants
result_labels_contam.add("contaminant")
result_labels_contam.add(CONTAM_LOW_COV)  # AA: not sure why the one above is literal and this one is named-constant.


def is_primary(s) -> bool:
    return s in result_labels_primary


def is_contam(s) -> bool:
    return s in result_labels_contam


def is_grey(s) -> bool:
    return not (is_contam(s) or is_primary(s))


class ResultGroup(Enum):
    CONTAM = 1
    PRIMARY = 2
    INCONCLUSIVE = -1


def result_group(s) -> ResultGroup:
    # 1 - contam, 2 - primary, -1 - inconclusive
    return (
             ResultGroup.CONTAM  if is_contam(s)
        else ResultGroup.PRIMARY if is_primary(s)
        else ResultGroup.INCONCLUSIVE
    )


#############################################################################
# Query-specific leq-id,len,masked-len,len-by-all, and up-to 4 taxa subrecords
@dataclass
class Record:
    """record"""

    seq_id      : str
    len         : int  # sequence length
    cvg_by_all  : int  # coverage length by all alignments
    tax_name    : str
    taxa        : List[Taxon]  # exactly (or up-to) 4?
    result      : str
    div         : str
    div_pct_cvg : int
    xtrachr_len : int  # raw extrachromoomal coverage (plastids, plasmids, mito)

    orig_id     : str
    start_pos   : int
    end_pos     : int
    whole       : bool
    gaps        : list
    # AA: Based on the code below, `gaps` is a 2-element list of bools
    # corresponding to whether the chimeric (~~) location abuts
    # the query start and stop respectively,
    # i.e. set to True iff there's NO gap at the corresponding end.
    # [True, True] if not chimeric loc.

    #########################################################################

    @staticmethod
    def id_to_range(seq_id: str):  # -> SeqRange

        SeqRange = namedtuple("SeqRange", "id, start, end, whole, gaps")

        tilda_count = seq_id.count("~")

        if tilda_count == 0:
            return SeqRange(seq_id, 0, 0, True, [True, True])

        if tilda_count in (1, 2):  # LR994621.1~46904608..46976173
            i1 = seq_id.split("~")
            i2 = i1[-1].split("..")
            return SeqRange(i1[0], int(i2[0]), int(i2[1]), False, [True, True])

        assert tilda_count == 3 # LR994621.1~20258968..20296541~~33292..35019

        i1 = seq_id.split("~")  # [LR994621.1, 20258968..20296541, , 33292..35019]
        i2 = i1[1].split("..")  # [20258968, 20296541]
        i3 = i1[3].split("..")  # [33292..35019]

        total_len = int(i2[1]) - int(i2[0]) + 1
        start_pos = int(i2[0]) + int(i3[0]) - 1
        end_pos   = int(i2[0]) + int(i3[1]) - 1

        eg = i3[0] == "1"
        lg = int(i3[1]) == total_len

        return SeqRange(i1[0], start_pos, end_pos, False, [eg, lg])


    #########################################################################
    @staticmethod
    def from_row(row: List[str]):
        assert len(row) >= 30 and row[4] == "|" and row[29] == "|"
        assert row[31] in result_labels_valid

        lens = [int(s) for s in row[2].split(",")]  # [transposons, low-complexity, conserved, n-runs, extrachromosomal]
        lens += [0] * (5 - len(lens))               # pad to len=5 (extrachromosomal is present only in newer outputs).
        assert len(lens) >= 5

        rec = Record(
            seq_id      = row[0],
            len         = int(row[1]),
            cvg_by_all  = int(row[3]),
            tax_name    = row[5],
            taxa        = [],
            result      = row[31],
            div         = row[32],
            div_pct_cvg = int(row[33]),
            xtrachr_len = lens[4],
            orig_id     = "",
            start_pos   = 0,
            end_pos     = 0,
            whole       = False,
            gaps        = [],
        )

        for i in [6, 12, 18, 24]:  # columns at which subrecords begin
            # Only consider taxa with scores at least some minimum noise-floor threshold (25)
            if row[i + 4] and int(row[i + 4]) > 25:
                rec.taxa.append(Taxon.from_row(row, i))

        if "~" in rec.seq_id:
            row_range     = Record.id_to_range(rec.seq_id)
            rec.orig_id   = row_range.id
            rec.start_pos = row_range.start
            rec.end_pos   = row_range.end
            rec.whole     = row_range.whole
            rec.gaps      = row_range.gaps
        else:
            rec.orig_id   = rec.seq_id
            rec.start_pos = 1
            rec.end_pos   = rec.len
            rec.whole     = True

        return rec


#############################################################################
# AA: looks like this should be a dataclass?
class Span:

    # AA: why do we have both __init__() and make_span()?

    def __init__(
        self,
        records: List[Record],
        result: str,
        orig_id: str,
        start_pos: int,
        end_pos: int,
        whole: bool
    ):
        self.records    = records
        self.result     = result
        self.orig_id    = orig_id
        self.start_pos  = start_pos
        self.end_pos    = end_pos
        self.whole      = whole

        self.rules: List[List[str]] = []

        self.total_len      = 0
        self.ungap_len      = 0
        self.result_len     = 0
        self.contam_len     = 0
        self.primary_len    = 0
        self.grey_len       = 0
        self.contam_len_lc  = 0
        self.result_cov     = 0
        self.top_tax_name   = ""
        self.max_div_name   = ""
        self.max_div_cov    = 0
        self.result_cov     = 0

    def __repr__(self):
        return f"Span(len(records)={len(self.records)}, result={self.result}, orig_id={self.orig_id}, start_pos={self.start_pos}, end_pos={self.end_pos}, top_tax_name={self.top_tax_name}, max_div_name={self.max_div_name})"

    #########################################################################
    @staticmethod
    def make_span(records: List[Record]):

        if len(records) == 0:
            return None

        span = Span(
            records,
            records[0].result,
            records[0].orig_id,
            records[0].start_pos,
            records[-1].end_pos,
            records[0].whole
        )

        span.rules = []

        span.recalc_span_stats()

        return span

    #########################################################################
    def merge_spans(self, others, force_contam=False):
        for span in others:
            self.records += span.records

        self.recalc_span_stats(force_contam)

    #########################################################################
    def recalc_span_stats(self, force_contam=False):
        seq_stat = SeqStat(self.orig_id)
        seq_stat.add_span(self, force_contam)
        seq_stat.finalize()

        self.result = seq_stat.result
        self.start_pos = seq_stat.start_pos
        self.end_pos = seq_stat.end_pos
        self.total_len = seq_stat.total_len
        self.ungap_len = seq_stat.ungap_len
        self.result_len = seq_stat.result_len
        self.contam_len = seq_stat.contam_len
        self.primary_len = seq_stat.primary_len
        self.grey_len = seq_stat.grey_len
        self.contam_len_lc = seq_stat.contam_len_lc
        self.result_cov = seq_stat.result_cov
        self.top_tax_name = seq_stat.top_tax_name
        self.max_div_name = seq_stat.max_div_name
        self.max_div_cov = seq_stat.max_div_cov

        return self


    #########################################################################
    @staticmethod
    def grey_side__GP_36740(left: Record, right: Record) -> Tuple[bool, str]:
        """
        This is putative re-implementation of grey_side(...) below,
        based on the following discussion in Slack:
        https://ncbi.slack.com/archives/C0277FK9FH7/p1697816519856089

        return (b:bool, s:str)
        where b = True if left is dominant based on raw-coverage(penalized-if-contam), else False
              s = "G01" if dominant is contam, else "G02"

        Note: definitions of G## from down below in this script:
        "G01": "group 'grey' rows with the high-coverage contaminant hits",
        "G02": "group 'grey' rows with the high coverage primary-div hits",
        "G03": "group 'grey' rows with the contaminant, as default",  # NB: unused in this version


        NB: Not sure if penalizing by 0.75 factor if contam is necessary
        in get_weight() below.
       
        I added that based on the following comment in the above thread:
        "or for a weak contam include the greys with the stronger primarys
        ... and the 75% checks are because if it was just 
        p.pident  > r.pident it would be fighting between a 49% and a 51%"
        
        So I interpret it as follows: For contam vs. primary case of 
        about the same raw coverage (I presume 'r.pident' was a mistype),
        we want primary to 'win' - i.e. "greys" to be grouped with it,
        hence we need to penalize contam by some factor less than 1).
        """

        assert is_contam(left.result) != is_contam(right.result)

        def get_weight(r):
            return r.div_pct_cvg * r.len * (0.75 if is_contam(r.result) else 1)

        left_is_dominant   = get_weight(left) > get_weight(right)
        dominant_is_contam = is_contam(left.result if left_is_dominant else right.result)

        return (left_is_dominant, "G01" if dominant_is_contam else "G02")


    # The original grey_side() implementation is below, which is diverted to the above
    # if enable_putative_fix is set. If the above implementation is approved,
    # then the one below can be replaced for good.

    #########################################################################
    # assuming records ordered
    #   p ... some greys ... r
    # returns: True if the greys should be grouped with p
    #          False if the greys should be grouped with r
    @staticmethod
    def grey_side(p: Record, r: Record, enable_putative_fix=int(os.getenv("GX_ACTION_REPORT_ENABLE_GP_36740", "0"))):

        if enable_putative_fix:
            return Span.grey_side__GP_36740(p, r)

        pending_goes_early = False
        rule = ""

        if is_contam(p.result):  # p is contam and r is primary
            if p.div_pct_cvg >= 75:
                pending_goes_early = False
                rule = "G01"

            elif r.div_pct_cvg >= 75:
                pending_goes_early = True
                rule = "G02"

            else:
                pending_goes_early = False
                rule = "G03"

        else:  # p is primary-div, r is contam
            if r.div_pct_cvg >= 75:
                pending_goes_early = True
                rule = "G01"

            elif p.div_pct_cvg >= 75:
                pending_goes_early = False
                rule = "G02"

            else:
                pending_goes_early = True
                rule = "G03"

        dprint(f"    p: {p.result} : {p.div_pct_cvg}")
        dprint(f"    r: {r.result} : {r.div_pct_cvg}")
        dprint(f"  pending goes early: {pending_goes_early}")
        dprint(f"   rule: {rule}")

        return pending_goes_early, rule

    #########################################################################
    @staticmethod
    def calc_spans(records: List[Record]):
        # make spans of similar records
        span_records = [[]]
        pending_greys = []
        rules = [[]]

        for r in records:
            dprint(f"span_records len:{len(span_records)}")
            dprint(f"last span len:{len(span_records[-1])}")
            dprint(f"pending_greys len:{len(pending_greys)}")
            dprint(f"curr result: {r.result}  start_pos:{r.start_pos}")
            dprint("----")

            if result_group(r.result) == ResultGroup.INCONCLUSIVE:
                pending_greys.append(r)
                continue

            if len(span_records[-1]) == 0:
                span_records[-1].extend(pending_greys)
                pending_greys.clear()
                span_records[-1].append(r)

            elif result_group(span_records[-1][-1].result) == result_group(r.result):
                span_records[-1].extend(pending_greys)
                pending_greys.clear()
                span_records[-1].append(r)

            else:  # transitioning to other state
                p = span_records[-1][-1]
                new_span = []
                pending_goes_early, rule = Span.grey_side(p, r)

                if pending_goes_early:
                    span_records[-1].extend(pending_greys)
                    rules[-1].append(rule)
                    rules.append([])
                else:
                    new_span.extend(pending_greys)
                    rules.append([rule])

                pending_greys.clear()
                new_span.append(r)
                span_records.append(new_span)

        span_records[-1].extend(pending_greys)
        pending_greys.clear()

        spans: List[Span] = []

        for si, s in enumerate(span_records):
            span = Span.make_span(s)

            if span is None:
                continue

            span.rules = rules[si]
            spans.append(span)

            dprint(
                f"{si+1}:\t{span.result}:\t{len(span.records)}\t{span.orig_id}:\t"
                f"{span.start_pos} .. {span.end_pos}\t{':'.join(span.rules)}  :  {span.result_len} after grey merge"
            )

        if len(spans) == 1:
            return spans

        mspans = []
        # i = 0
        # j = len(spans) - 1

        Range = namedtuple("Range", "i j")  # AA: based on the code below, j seems to be inclusive.
        poss_ranges = []

        # seed poss_ranges with ranges that exist between primary_len spans of >50K
        big_primary_spans = [i for (i, s) in enumerate(spans) if s.primary_len >= 50000]

        if len(big_primary_spans) > 0:
            pc = -1

            for bps in big_primary_spans:

                if bps <= pc + 1:
                    pc = bps
                    continue

                poss_ranges.append(Range(pc + 1, bps - 1))
                pc = bps

            if pc < len(spans):
                poss_ranges.append(Range(pc, len(spans) - 1))

            if poss_ranges[0].i > 0:
                poss_ranges.insert(0, Range(0, poss_ranges[0].i - 1))
        else:
            poss_ranges.append(Range(0, len(spans) - 1))

        # bps itself can be missing from poss_ranges, downstream incorrect merge if contam regions on both sides.
        # see if it is included, add if necesary without breaking existing logic above. eg GCF_010730955.1
        for bps in big_primary_spans:
            found = False

            for rr in poss_ranges:
                if bps in rr:
                    found = True
                    break

            if not found:
                poss_ranges.append(Range(bps, bps))

        poss_ranges.sort(key=lambda r: r.i)

        dprint(f"{len(spans)}")
        dprint(f"bps: {big_primary_spans} ")
        dprint(f"prs: {poss_ranges}")

        while poss_ranges:
            p = poss_ranges.pop()
            dprint(f"merge loop : {p.i} , {p.j}")

            if p.i > p.j:
                continue

            if p.i == p.j:
                dprint(f" Just keeping span: {p.i} \t {spans[p.i].start_pos} .. {spans[p.i].end_pos}")
                mspans.append(spans[p.i])
                continue

            i = p.i
            j = p.j

            while i < j:
                rr = list(itertools.chain.from_iterable(x.records for x in spans[i : j + 1]))
                # c = i
                # rr = []
                # while c <= j:
                #    if global_debug:
                #        s = spans[c]
                #        eprint(
                #            f"In test merge: ({i}, {j}) \t {s.start_pos}..{s.end_pos} \t "
                #            f"{s.contam_len}:{s.primary_len} \t {s.result}"
                #        )
                #    rr.extend(spans[c].records)
                #    c = c + 1

                rrs = Span.make_span(rr)
                if (
                    rrs.contam_len > 0
                    and rrs.primary_len > 0
                    and int(rrs.contam_len * 0.30) > rrs.primary_len
                    and rrs.primary_len < 50000
                ):  # pylint: disable=C0301
                    dprint(f"Merging: {i}, {j} \t " f"  lens: {rrs.contam_len} : {rrs.primary_len}")
                    mspans.append(rrs)
                    poss_ranges.append(Range(p.i, i - 1))
                    poss_ranges.append(Range(j + 1, p.j))
                    break

                if spans[i].result_len >= spans[j].result_len:
                    i += 1
                else:
                    j -= 1


            if i == j:
                dprint(f" Just keeping spans: {p.i} , {p.j} \t {spans[p.i].start_pos} .. {spans[p.j].end_pos}")
                mspans.extend(spans[p.i : p.j + 1])

                for si, mspan in enumerate(mspans):
                    dprint(
                        f"{si+1}:\t{mspan.result}:\t{len(mspan.records)}\t{mspan.orig_id}:\t"
                        f"{mspan.start_pos} .. {mspan.end_pos}\t{':'.join(mspan.rules)}:\t"
                        f"{mspan.contam_len}:{mspan.primary_len} i {i} j {j} after merging"
                    )

        mspans.sort(key=lambda x: x.start_pos)
        spans = mspans

        for si, span in enumerate(spans):
            dprint(
                f"{si+1}:\t{span.result}:\t{len(span.records)}\t{span.orig_id}:\t"
                f"{span.start_pos} .. {span.end_pos}\t{':'.join(span.rules)}:\t"
                f"{span.contam_len}:{span.primary_len}"
            )
        return spans


#############################################################################
class SeqStat:
    """Class for statistics on sequence or span"""
    def __init__(self, seq_id):
        self.seq_id = seq_id
        self.start_pos = -1
        self.end_pos = -1
        self.total_len = 0

        self.ungap_len = 0
        self.result_len = 0
        self.contam_len = 0
        self.primary_len = 0
        self.grey_len = 0
        self.contam_len_lc = 0

        self.result = ""
        self.any_result = ""

        self.top_tax_name = "n/a"

        self.max_div_name = ""
        self.max_div_cov = 0

        self.result_cov = 0

        self.taxname_lens: Dict[str, int] = defaultdict(int)

        self.div_lens: Dict[str, int] = defaultdict(int)

    def __repr__(self):
        return f"SeqStat(id={self.seq_id},result={self.result}, start_pos={self.start_pos}, end_pos={self.end_pos}, top_tax_name={self.top_tax_name}, max_div_name={self.max_div_name}, max_div_cov={self.max_div_cov}, taxname_lens={self.taxname_lens})"

    #########################################################################
    def add_span(self, span: Span, force_contam=False):
        dsprint(self.seq_id, f"adding span of {len(span.records)} records")

        for r in span.records:
            if result_group(r.result) != ResultGroup.INCONCLUSIVE:
                self.result = r.result
                break

        # covers when the rows are all only grey
        if self.result == "":
            # Save first encountered result till finalize
            if self.any_result == "":
                self.any_result = span.records[0].result

        if self.start_pos == -1:
            self.start_pos = span.records[0].start_pos
        else:
            self.start_pos = min(self.start_pos, span.records[0].start_pos)

        self.end_pos = max(self.end_pos, span.records[-1].end_pos)
        self.total_len = self.end_pos - self.start_pos + 1

        for r in span.records:
            self.ungap_len = self.ungap_len + r.len
            r_len = r.len * r.div_pct_cvg // 100
            rg = result_group(r.result)

            if rg == ResultGroup.CONTAM:
                self.contam_len = self.contam_len + r_len
            elif rg == ResultGroup.PRIMARY:
                self.primary_len = self.primary_len + r_len
            elif rg == ResultGroup.INCONCLUSIVE:
                self.grey_len = self.grey_len + r_len

        if (force_contam or self.contam_len >= self.primary_len) and self.contam_len > 0:
            self.result_len = self.contam_len
            self.result = "contaminant"

        elif self.primary_len > self.contam_len:
            self.result_len = self.primary_len
            self.result = "primary-div"

        else:
            self.result_len = self.grey_len
            self.result = "inconclusive"  # representative of any greys

        self.result_cov = int((self.result_len / self.ungap_len) * 100)

        taxname_lens = self.taxname_lens
        for r in span.records:
            if r.taxa and (is_contam(r.result) or (is_primary(r.result) and not force_contam)):
                taxname_lens[r.tax_name] += r.taxa[0].cvg_by_tax

        div_lens = self.div_lens

        for r in span.records:
            # if result_group(r.result) != ResultGroup.INCONCLUSIVE:
            if is_contam(r.result) or (is_primary(r.result) and not force_contam):
                div_lens[r.div] += r.len * r.div_pct_cvg // 100

        if div_lens:
            max_div_name = max(div_lens, key=div_lens.get)

            for r in span.records:
                if r.result == "low-coverage" and r.div == max_div_name:
                    div_lens[r.div] += r.len * r.div_pct_cvg // 100

    #########################################################################
    def finalize(self):
        if not self.result:
            self.result = self.any_result

        if self.taxname_lens:
            self.top_tax_name = max(self.taxname_lens, key=self.taxname_lens.get)

        div_lens = self.div_lens

        if div_lens:
            self.max_div_name = max(div_lens, key=div_lens.get)

            if self.ungap_len > 0:
                self.max_div_cov = int((div_lens[self.max_div_name] / float(self.ungap_len)) * 100)

        dsprint(self.seq_id, "stat is", self)

#############################################################################
@dataclass
class Action:
    seq_id       : str
    start_pos    : int
    end_pos      : int
    len          : int
    action       : str
    div          : str
    agg_cov      : int
    top_tax_name : str

    span         : Optional[Span]  # AA: can be None (see FIXME/"A00")
    rule         : str

    #########################################################################
    def tsv(self, print_rules) -> str:
        return "\t".join(
            [
                self.seq_id.replace("{TILDE}", "~"),  # GP-35596
                str(self.start_pos),
                str(self.end_pos),
                str(self.len),
                self.action,
                self.div,
                str(self.agg_cov),
                self.top_tax_name,
            ] + ([":".join([*self.span.rules, self.rule])] if print_rules and self.span else [])
        )


@dataclass
class ExcludeThresholds:
    good_pct: int = 30
    good_len: int = 50000
    good_abs_len: int = 200


re_endosymbiont = re.compile(r"wolbachia", re.I)

def keep(rule="???"):
    return ("KEEP", rule) if global_debug else ("", "")


@dataclass
class Sequence:
    """sequence"""

    seq_id: str
    seq_len: int
    ungap_len: int
    spans: List[Span]
    actions: list

    #########################################################################
    @staticmethod
    def calc_sequence(records):

        if len(records) == 0:
            return None

        seq_id = records[0].orig_id
        dsprint(seq_id, f"calc_sequence for {seq_id}")
        seq_len = 0
        ungap_len = 0

        actions = []
        spans = Span.calc_spans(records)

        # row-by-row data
        for r in records:
            seq_len = max(r.end_pos, seq_len)
            ungap_len = ungap_len + r.len

        seq = Sequence(seq_id, seq_len, ungap_len, spans, actions)
        dsprint(seq_id, "calc_sequence result", seq)
        return seq

    #########################################################################
    def discard_small_primary_islands(self, thr: ExcludeThresholds):
        """ Expand trims as long as primary divs are smaller than percentage of contamination """

        if not global_discard_small_primaries:
            return

        spans = self.spans
        seq_id = self.seq_id
        trim_index_front = Sequence.merge_trim_index(spans, thr)
        trim_index_back = len(spans) - Sequence.merge_trim_index(reversed(spans), thr)

        if trim_index_front > 0:
            dsprint(seq_id,
                f"Seqid {seq_id} front trim from {spans[0].start_pos} "
                f"to {spans[trim_index_front].end_pos}"
            )
            # Index is inclusive, so to get range [i, j] we need to use Python expression a[i:j+1]
            trim_index_front += 1

        if trim_index_back < len(spans):
            dsprint(seq_id,
                f"Seqid {seq_id} back trim from {spans[trim_index_back-1].start_pos} "
                f"to {spans[-1].end_pos}"
            )
            # Same here
            trim_index_back -= 1

        trim_front = spans[0:trim_index_front]

        if trim_front:
            merged_span = trim_front[0]
            merged_span.merge_spans(trim_front[1:], global_top_contam_organism)
            merged_span.result = "contaminant"
            trim_front = [merged_span]

        trim_back = spans[trim_index_back : len(spans)]

        if trim_back:
            merged_span = trim_back[0]
            merged_span.merge_spans(trim_back[1:], global_top_contam_organism)
            merged_span.result = "contaminant"
            trim_back = [merged_span]

        self.spans = trim_front + spans[trim_index_front:trim_index_back] + trim_back

    #########################################################################
    def range_has_gaps(self, left, right):

        for j in range(left, right):
            first = self.spans[j]
            second = self.spans[j+1]

            if first.end_pos < second.start_pos - 1:
                if (
                       (first.result  == "primary-div" and is_contam(second.result))
                    or (second.result == "primary-div" and is_contam(first.result))
                ):
                    return True

        return False

    #########################################################################
    def find_integrants(self):
        """Find contamination elements which are most probably integrants
           Following is believed to be an integrant:
           Euk virus on euks and prok virus on prok, prok on euk.
           Integrant must be flanked by compatible primary divs without gaps
           between primaries and integrant. Gaps between inconclusive and primary/integrant
           are allowed. Also, gaps are not checked for endosymbionts Wolbachia and Rickettsia.
           Returns:
               list of pairs of integrant spans, primary divs included
        """
        spans = self.spans
        primary_div_indices = [n for n in range(len(spans)) if spans[n].result == "primary-div"]

        if len(primary_div_indices) < 2:
            return []

        # eprint(f"{self.seq_id}, primary divs: {', '.join(map(str, primary_div_indices))}")
        integrants = []

        for i in range(len(primary_div_indices)-1):
            left_flank = primary_div_indices[i]
            right_flank = primary_div_indices[i+1]

            if right_flank - left_flank < 2:
                continue

            left_flank_div = spans[left_flank].max_div_name
            right_flank_div = spans[right_flank].max_div_name

            if is_prok(left_flank_div) != is_prok(right_flank_div):
                continue

            # Check for compatibility: euk virus on euks, prok virus on proks
            # prok on euk endosymbiont (currently checks specifically for Wolbachia and Rickettsia)
            compat = True
            for j in range(left_flank+1, right_flank):
                s = spans[j]

                if s.result == "inconclusive":
                    continue

                if re_endosymbiont.search(s.top_tax_name) and is_euk(left_flank_div):
                    # Valid endosymbiont
                    # eprint(f"{self.seq_id} [{left_flank+1}:{right_flank-1}] div {s.max_div_name} tax {s.top_tax_name} Valid endosymbiont")
                    continue

                if is_virus(s.max_div_name):
                    if self.range_has_gaps(left_flank, right_flank) or not can_host_virus(left_flank_div, s.max_div_name):
                        # eprint("Virus is of wrong kind")
                        compat = False
                        break

                else:
                    # Neither a virus nor an endosymbiont
                    # eprint(f"{self.seq_id} [{left_flank+1}:{right_flank-1}] div {s.max_div_name} tax {s.top_tax_name} Contamination is of wrong kind")
                    compat = False
                    break

            if not compat:
                continue

            # eprint(f"Integrant candidate: {self.seq_id} [{left_flank+1}:{right_flank-1}] {','.join(map(lambda x: x.max_div_name, spans[left_flank+1:right_flank]))}")
            integrants.append((left_flank, right_flank))

        return integrants

    #########################################################################
    def check_for_endosymbiotic_chimeras(self, primary_div_name: str) -> bool:
        # To reduce number of false positives we don't report bacteria-in-insect-or-nematode chimeras if
        # bacteria is Wolbachia or Rickettsia. Chimera is defined here very broadly - any combination
        # of insect or nematode primary div, repeat, low-coverage; and Wol/Ric contaminant
        if len(self.spans[0].records) == 0:
            return False

        seq_id = self.spans[0].records[0].seq_id
        dsprint(seq_id, f"Checking for endosymbiont chimera {seq_id} with ptimary div {primary_div_name}")

        if primary_div_name not in {"anml:insects", "anml:nematodes"}:
            return False

        is_chimeric = True
        host_present = False
        endosymbiont_present = False

        for span in self.spans:
            for r in span.records:
                dsprint(seq_id, f"Record {r.seq_id}, tax_name={r.tax_name}, result={r.result}, div={r.div}")

                if re_endosymbiont.search(r.tax_name) and is_contam(r.result):
                    # endosymbiotic contaminant
                    endosymbiont_present = True
                    continue

                if (is_primary(r.result) or r.result in {"repeat", "low-coverage"}) and r.div in {"anml:insects", "anml:nematodes"}:
                    # insect or nematode host
                    host_present = True
                    continue

                is_chimeric = False
                break

        dsprint(seq_id, f"is_chimeric={is_chimeric}, host_present={host_present}, endosymbiont_present={endosymbiont_present}")
        return is_chimeric and host_present and endosymbiont_present

    #########################################################################
    def calc_actions(
        self,
        primary_div_name: str,
        low_cov_mode: bool,
        pa_same_kingdom_contam_demote: bool,
        is_same_kingdom_demote: bool,
        thr: ExcludeThresholds
    ):
        spans     = self.spans
        actions   = self.actions
        seq_id    = self.seq_id
        seq_len   = self.seq_len
        ungap_len = self.ungap_len


        dsprint(seq_id, f"good_len={thr.good_len}, good_pct={thr.good_pct}, good_abs_len={thr.good_abs_len}")

        if len(spans) == 0:
            actions.append(Action(seq_id, 0, 0, 0, "FIXME", "", 0, "", None, "A00"))

        elif len(spans) == 1:
            s = spans[0]
            action_str = ""
            act_start = 1
            act_end = seq_len

            dsprint(seq_id, f"825: id: {self.seq_id},  low_cov_mode: {low_cov_mode}, isk: {is_same_kingdom(primary_div_name, s.max_div_name)}, pdn: {primary_div_name}, mdn: {s.max_div_name}, mdc: {s.max_div_cov}, result={s.result}")


            if (  # GP-36649
                is_primary(s.result)
                and is_euk(primary_div_name)
                and sum(r.xtrachr_len for span in spans for r in span.records) > 0.8 * seq_len
            ):
                action_str = "ORGANELLE"
                rule = "S04"  # Extrachromosomal coverage > 80% (plastid or mito)

            elif is_primary(s.result):
                action_str = "ACCEPT"
                rule = "S02"

            elif s.max_div_cov <= global_low_coverage_threshold:
                # Should we ACCEPT low coverage records or MEH?
                action_str = "ACCEPT"
                rule = "coverage too low to decide"

            elif is_contam(s.result):
                # elif is_euk(primary_div_name) and re_endosymbiont.search(s.top_tax_name):
                #     action_str = "INFO"
                #     rule = "Endosymbiont"
                # GP-34699 - move logic for prok virus in prok upstream
                # if is_prok(primary_div_name) and is_prok_virus(s.max_div_name):
                #     action_str = "ACCEPT"
                #     rule = "prok virus in prok"
                # elif ...

                if pa_same_kingdom_contam_demote and is_prok_or_arch(s.max_div_name):
                    # pa_same_kingdom_contam_demote implies is_prok_or_arch(primary_div_name)
                    assert is_prok_or_arch(primary_div_name)
                    action_str = "REVIEW_RARE"
                    rule = f"p/a same kingdom {global_pa_same_kingdom_threshold}% treshold"

                elif is_same_kingdom_demote and is_same_kingdom(s.max_div_name, primary_div_name):
                    action_str = "ACCEPT"
                    rule = "same kingdom contam ignored for at least TSA inputs"

                elif self.check_for_endosymbiotic_chimeras(primary_div_name):
                    action_str = "INFO"
                    rule = "endosymbiotic chimera"
                    dsprint(seq_id, f"Endosymbiotic chimera {seq_id}")

                elif is_mostly_review_contam(s.records):
                    action_str = "REVIEW"
                    rule = "low-coverage review"
                    dsprint(seq_id, f"Low-coverage Review {seq_id}")

                else:
                    action_str = "EXCLUDE"
                    rule = "S01"
                    dsprint(seq_id, f"852: {action_str} : isk: {is_same_kingdom(primary_div_name, s.max_div_name)}, tl: {s.total_len}, cov: {s.max_div_cov}, mdn: {s.max_div_name}")

                    if is_same_kingdom(primary_div_name, s.max_div_name) and s.max_div_cov < 20:
                        action_str = "REVIEW"
            else:
                action_str = "MEH"
                rule = "S03"

            if global_debug or action_str in ("EXCLUDE", "REVIEW", "REVIEW_RARE", "ORGANELLE", "INFO"):
                if action_str in ("EXCLUDE", "REVIEW", "REVIEW_RARE", "ORGANELLE"):
                    s.recalc_span_stats(global_top_contam_organism)

                actions.append(
                    Action(
                        seq_id,
                        act_start,
                        act_end,
                        seq_len,
                        action_str,
                        s.max_div_name,
                        s.max_div_cov,
                        s.top_tax_name,
                        s,
                        rule,
                    )
                )

        elif len(spans) >= 2:
            # mixed, but the good parts are bad enough to just exclude it anyway
            keep_len = 0  # ungapped total length of good spans
            keep_result_len = 0  # total coverage by primary div of good spans
            first_contam_index = len(spans)
            i = 0

            for s in spans:

                if is_primary(s.result):
                    keep_len = keep_len + s.ungap_len
                    keep_result_len = keep_result_len + s.result_len

                elif is_contam(s.result):
                    first_contam_index = min(first_contam_index, i)

                i = i + 1

            dsprint(seq_id, f"keep_len={keep_len}, keep_result_len={keep_result_len}, ungap_len={ungap_len}")

            # GP-32655 if the total primary div length is too small, discard the whole sequence
            # Total primary div length, keep_len can be either too small in absolute sense (thr.good_abs_len)
            # or if it is smaller than some larger threshold (thr.good_len) and smaller than some
            # percentage (thr.good_pct) of total ungapped length for the sequence (ungap_len)
            if (keep_len < thr.good_abs_len) or (
                first_contam_index < len(spans)
                and (keep_len < thr.good_len)
                and (keep_len / float(ungap_len) * 100) < thr.good_pct
            ):
                dsprint(seq_id, f"Primary divs are too short for {seq_id}, classify as contamination")

                if first_contam_index < len(spans):
                    s = spans[first_contam_index]

                action_str = ""
                seq_stat = SeqStat(seq_id)

                # GP-33831 recalculate contamination coverage (or all coverage?)
                for s in spans:
                    seq_stat.add_span(s, global_top_contam_organism)

                seq_stat.finalize()
                pct_cvg = seq_stat.max_div_cov
                max_div_name = seq_stat.max_div_name
                top_tax_name = seq_stat.top_tax_name
                dsprint(seq_id, f"910: id: {self.seq_id}, as: {action_str}, tl: {s.total_len}, sr: {s.result}, cov: {s.max_div_cov}, low_cov_mode: {low_cov_mode}, "
                        f"isk: {is_same_kingdom(primary_div_name, s.max_div_name)}, pdn: {primary_div_name}, mdn: {s.max_div_name}, mdc: {s.max_div_cov}")
                # For very low coverage organism we don't have enough information to exclude even bad
                # spans, so we accept them

                if pct_cvg <= global_low_coverage_threshold:
                    action_rule = "M05"

                    if global_debug:
                        action_str = "ACCEPT" # Should we use MEH here?
                    else:
                        action_str = ""

                elif pa_same_kingdom_contam_demote and is_prok_or_arch(max_div_name):
                    # pa_same_kingdom_contam_demote implies is_prok_or_arch(primary_div_name)
                    assert is_prok_or_arch(primary_div_name)
                    action_str = "REVIEW_RARE"
                    action_rule = f"p/a same kingdom {global_pa_same_kingdom_threshold}% treshold"

                elif s.result == CONTAM_LOW_COV:
                    action_str = "REVIEW"
                    action_rule = "low-coverage review"
                    dsprint(seq_id, f"Low-coverage Review {seq_id}")

                else:
                    action_str = "EXCLUDE"
                    action_rule = "M01"
                    #GP-35454: prok/archea treated as same kingdom here

                    if (is_same_kingdom_pa(primary_div_name, seq_stat.max_div_name) and seq_stat.max_div_cov < 20):
                        action_str = "REVIEW"
                        action_rule = "same kingdom contaminant with coverage less than 20%"

                if action_str:
                    actions.append(
                        Action(
                            seq_id,
                            1,
                            seq_len,
                            seq_len,
                            action_str,
                            max_div_name,
                            pct_cvg,
                            top_tax_name,
                            s,
                            action_rule,
                        )
                    )
            else:
                # actual mixed, span by span
                dsprint(seq_id, f"Mixed sequence {seq_id}")
                self.discard_small_primary_islands(thr)
                # Use actual spans, not cached old variable
                spans = self.spans

                integrant_set = set()
                if global_find_integrants:

                    integrants = self.find_integrants()

                    for b, e in integrants:
                        for i in range(b + 1, e):
                            integrant_set.add(i)

                for i, s in enumerate(spans):
                    trim_action = ""
                    trim_rule = ""
                    dprint(f"1077: mdn {s.max_div_name}, mdc: {s.max_div_cov}, start {s.start_pos}, end {s.end_pos}")

                    if is_primary(s.result):
                        # Keep spans which match main org or have coverage
                        # <= 2% - we don't have enough info to reject them
                        # GP-33724
                        # Should we KEEP it, or report as MEH in debug mode?
                        trim_action, trim_rule = keep()

                    elif s.max_div_cov <= global_low_coverage_threshold:
                        trim_action, trim_rule = keep("coverage too low to decide")

                    elif is_same_kingdom_demote and is_same_kingdom(s.max_div_name, primary_div_name):
                        trim_action, trim_rule = keep("same kingdom contam ignored for at least TSA inputs")

                    elif is_contam(s.result):

                        if global_debug_seq_id == seq_id:
                            eprint(f"965: id:{self.seq_id} , tl {s.total_len} : cov {s.max_div_cov} ,"
                                   f" start {s.start_pos}, end {s.end_pos},"
                                   f" low_cov_mode: {low_cov_mode} , isk: {is_same_kingdom(primary_div_name, s.max_div_name)} ,"
                                   f" pdn: {primary_div_name} , mdn: {s.max_div_name}")

                        trim_action = "FIX"
                        trim_rule = "M03"

                        if s.start_pos == 1 or s.end_pos == seq_len:
                            trim_action = "TRIM"
                            trim_rule = "M02"

                        # GP-34699 - move logic for prok virus in prok upstream
                        # if is_prok(primary_div_name) and is_prok_virus(s.max_div_name):
                        #     trim_action, trim_rule = keep("prok virus in prok")
                        # elif ...
                        if is_euk(primary_div_name) and re_endosymbiont.search(s.top_tax_name):
                            trim_action = "INFO"
                            trim_rule = "Endosymbiont"

                        elif i in integrant_set:
                            if is_prok(s.max_div_name) or re_endosymbiont.search(s.top_tax_name):
                                trim_action = "INFO"
                                trim_rule = "Endosymbiont"

                            elif is_virus(s.max_div_name):
                                trim_action, trim_rule = keep("integrant virus")

                        elif is_mostly_review_contam(s.records):
                            trim_action = "REVIEW"
                            trim_rule = "low-coverage review"

                        else:
                            assert trim_action in {"FIX", "TRIM"}
                            #GP-35454 treat prok and arch as same kingdom here

                            if is_same_kingdom_pa(primary_div_name, s.max_div_name):
                                if s.total_len * s.max_div_cov / 100 < 10000:
                                    trim_action, trim_rule = keep("same kingdom contaminant with effective length less than 10000")

                                elif pa_same_kingdom_contam_demote and is_prok_or_arch(s.max_div_name):
                                    # pa_same_kingdom_contam_demote implies is_prok_or_arch(primary_div_name)
                                    assert is_prok_or_arch(primary_div_name)
                                    trim_action = "REVIEW_RARE"
                                    trim_rule = f"p/a same kingdom {global_pa_same_kingdom_threshold}% treshold"

                                else:
                                    trim_action = "REVIEW"
                                    trim_rule = "same kingdom contaminant with effective length more than 10000 bases"

                            elif pa_same_kingdom_contam_demote and is_prok_or_arch(s.max_div_name):
                                # pa_same_kingdom_contam_demote implies is_prok_or_arch(primary_div_name)
                                assert is_prok_or_arch(primary_div_name)
                                trim_action = "REVIEW_RARE"
                                trim_rule = f"p/a same kingdom {global_pa_same_kingdom_threshold}% treshold"

                        if trim_action == "REVIEW" and low_cov_mode and is_same_kingdom(primary_div_name, s.max_div_name):
                            # Don't review in low coverage mode if the div is in the same kingdom as primary div
                            trim_action, trim_rule = keep("same kingdom low coverage contaminant")

                    if trim_action:
                        actions.append(
                            Action(
                                seq_id,
                                s.start_pos,
                                s.end_pos,
                                seq_len,
                                trim_action,
                                s.max_div_name,
                                s.max_div_cov,
                                s.top_tax_name,
                                s,
                                trim_rule,
                            )
                        )

        seq = Sequence(seq_id, seq_len, ungap_len, spans, actions)
        dsprint(seq_id, "sequence with actions added", seq)
        return seq

    #########################################################################
    @staticmethod
    def merge_trim_index(spans: List[Span], thr: ExcludeThresholds) -> int:
        """Find index of span which forms the merged trim according to
        criteria of GP-32594:
        if primary div is very short (ungapped len < good_abs_len) or
        short (ungapped len < good_len) and
        accumulated ungapped length of primary divs is shorter than good_pct of
        accumulated ungapped length of fix regions then we merge it over to a
        single trim region.
        """
        trim_index = 0
        trim_fix_ungap_len = 0
        trim_div_ungap_len = 0
        s_prev = None

        for i, s in enumerate(spans):

            if is_contam(s.result):
                trim_fix_ungap_len += s.ungap_len

                if (
                    i > 0
                    and (not s_prev or s_prev.ungap_len >= thr.good_abs_len)
                    and (trim_div_ungap_len / float(trim_fix_ungap_len) * 100 >= thr.good_pct)
                ):
                    break

                s_prev = None
                trim_index = i

            elif is_primary(s.result):

                if s.ungap_len >= thr.good_len:
                    break

                s_prev = s
                trim_div_ungap_len += s.ungap_len

        return trim_index

    #########################################################################
    @staticmethod
    def update_low_covs(seq: "Sequence", lc_review_divs: List) -> "Sequence":
        new_records = []

        for span in seq.spans:
            new_records.extend(span.records)

        # GP-34332, disable low-cov promotions, for now
        # GP-34277, add back, but with a list of divs instead of one
        for r in new_records:
            if r.result == "low-coverage" and r.div in lc_review_divs:
                #####and not is_same_kingdom( r.div):
                if r.div_pct_cvg >= 10:
                    r.result = CONTAM_LOW_COV

        return Sequence.calc_sequence(new_records)


#############################################################################
def calc_div_stats(seqs: List, primary_div_name: str):

    # AA: this can be a single dict: div -> dataclass
    div_lens                = defaultdict(int)
    div_contam_lens         = defaultdict(int)
    div_scaled_lens         = defaultdict(int)
    lc_div_lens             = defaultdict(int)
    div_contam_seq_counts   = defaultdict(int)
    div_seq_counts          = defaultdict(int)
    div_rec_counts          = defaultdict(int)
    div_cov_sums            = defaultdict(int)

    contam_divs         = set()
    lc_review_divs      = set()
    primary_div_len     = 0
    ungap_len           = 0
    total_gap_len       = 0
    primary_div_count   = 0
    pa_contam_lens_sum  = 0
    pa_contam_count_sum = 0

    for seq in seqs:
        seq_contam_divs = set()
        ungap_len = ungap_len + seq.ungap_len
        seq_divs = set()

        for span in seq.spans:
            for r in span.records:
                total_gap_len += r.len

                if is_contam(r.result) :
                    if r.orig_id == r.seq_id:  # not chimeric
                        contam_divs.add(r.div)
                        div_scaled_lens[r.div] += r.len * r.div_pct_cvg // 100
                        div_contam_lens[r.div] += r.len

                        div_cov_sums[r.div] += r.div_pct_cvg
                        div_rec_counts[r.div] += 1

                        dprint(f"div: {r.div} , rslt: {r.result}, dsl: {div_scaled_lens[r.div]}, dcl: {div_contam_lens[r.div]}, dl_c: {(div_scaled_lens[r.div]/div_contam_lens[r.div])}  || dcs: {div_cov_sums[r.div]}, drc: {div_rec_counts[r.div]}, dc%: {div_cov_sums[r.div]/div_rec_counts[r.div]}  ")
                        seq_contam_divs.add(r.div)

                if result_group(r.result) != ResultGroup.INCONCLUSIVE:
                    div_lens[r.div] += r.len
                    seq_divs.add(r.div)

                if result_group(r.result) != ResultGroup.INCONCLUSIVE or r.result == "low-coverage":
                    #lc_div_lens[r.div] += r.len * r.div_pct_cvg // 100
                    lc_div_lens[r.div] += r.len

                if is_contam(r.result) and is_prok_or_arch(r.div):
                    pa_contam_lens_sum += r.len
                    pa_contam_count_sum += 1

        for d in seq_divs:
            div_seq_counts[d] += 1

        for d in seq_contam_divs:
            div_contam_seq_counts[d] += 1


    # If primary div is not present in the header of the input file
    # we use primary div with most coverage
    if not primary_div_name:
        primary_div_name = max(div_lens, key=div_lens.get)

    if primary_div_name != "" and primary_div_name in div_lens:
        primary_div_len = div_lens[primary_div_name]
        primary_div_count = div_seq_counts[primary_div_name]
        div_lens.pop(primary_div_name)
        div_seq_counts.pop(primary_div_name)


    #  compute sum of contaminant calls for each division (# of sequences and total length).
    #    Divisions with >= 250 kb total row lengths are candidates from a given division for adding more REVIEW
    dprint(primary_div_name, primary_div_len, primary_div_count)
    for d in contam_divs:
        #div_avg_cov = div_cov_sums[d]/div_rec_counts[d]
        #dprint("contam: ", d, ", dcl: ", div_contam_lens[d], ", drc:", div_rec_counts[d], " dsc:", div_seq_counts[d], " dcsc:", div_contam_seq_counts[d],
        #       ", dcs:", div_cov_sums[d], ", dac:", div_avg_cov, "  alt: ", ungap_len, div_scaled_lens[d], (div_scaled_lens[d]/ungap_len) )
        if div_contam_lens[d] >= 250000 and not is_same_kingdom(primary_div_name, d):
            lc_review_divs.add(d)

    pa_same_kingdom_contam_demote = False
    if is_prok_or_arch(primary_div_name):
        if int((pa_contam_lens_sum / float(total_gap_len))*100) < global_pa_same_kingdom_threshold:
            pa_same_kingdom_contam_demote = True

    primary_div_cov = int((primary_div_len / float(ungap_len)) * 100)

    return {
        "ungap_len": ungap_len,
        "primary_div_name": primary_div_name,
        "primary_div_len": primary_div_len,
        "primary_div_cov": primary_div_cov,
        "primary_div_count": primary_div_count,
        "pa_same_kingdom_contam_demote": pa_same_kingdom_contam_demote,
        "lc_review_divs": lc_review_divs
    }


#############################################################################
def fix_action(action: Action, prefix: str) -> Action:
    action_str = "ACCEPT" if action.start_pos == 1 and action.end_pos == action.len else "KEEP"
    return Action(
        action.seq_id,
        action.start_pos,
        action.end_pos,
        action.len,
        action_str,
        action.div,
        action.agg_cov,
        action.top_tax_name,
        action.span,
        prefix + action.rule
    )


#############################################################################
def no_review_without_excludes_rule(seqs: List):
    # JIRA: GP-34021
    # if a div does not have any EXCLUDE rows
    # then there is not enough certainty about the divs
    # to make a REVIEW call either.
    # We can't just drop the action - we need to convert it to ACCEPT/KEEP for debug purposes

    divs_with_exclude = set((
        a.div
        for seq in seqs for a in seq.actions
        if a.action == "EXCLUDE"
    ))

    if global_debug:
        for seq in seqs:
            seq.actions = [
                (
                    a if a.action != "REVIEW" or a.div in divs_with_exclude
                    else fix_action(a, "No EXCLUDE for div, ")
                )
                for a in seq.actions
            ]
    else:
        for seq in seqs:
            seq.actions = [
                a
                for a in seq.actions
                if a.action != "REVIEW" or a.div in divs_with_exclude
            ]


#############################################################################
def production_review(seqs: List, ungap_len: int, build_name: str):
    # JIRA: GP-34563
    # contamination above threshold generates ancillary data for in-house review
    # thresholds: 1. FIX or TRIM 2. total contam length 3. contam division count 4. contam kingdom count
    CONTAM_KINGDOM_THRESH = 1
    CONTAM_DIV_THRESH = 4
    FT_LEN_THRESH = 50000
    ERFT_LEN_THRESH = 50000000
    CONTAM_LEN_THRESH = 0.2

    ft_len = []
    erft_len =[]
    contam_divs = set()
    contam_ks = set()
    for seq in seqs:
        for a in seq.actions:
            delta = a.end_pos -a.start_pos

            if a.action in ["FIX","TRIM"]:
                ft_len.append(delta)

            if a.action in ["EXCLUDE","REVIEW"]:
                erft_len.append(delta)
                contam_divs.add(a.div)

    for d in contam_divs:
        contam_ks.add(d[:5])

    err_string = ""
    if len(contam_ks) > CONTAM_KINGDOM_THRESH :
        err_string += "contaminant kingdoms; "

    if len(contam_divs) > CONTAM_DIV_THRESH :
        err_string += "contaminant divisions; "

    if sum(ft_len) > FT_LEN_THRESH or len(ft_len) > 2 :
        err_string += "fix+trim lengths; "

    if sum(erft_len) > ERFT_LEN_THRESH or sum(erft_len)/ungap_len > CONTAM_LEN_THRESH :
        err_string += "exclude+review lengths; "

    prod_warning = {}

    if err_string != "":
        prod_warning["build_name"] = build_name
        prod_warning["GX taxonomy analysis warning"] = err_string

    return prod_warning


#############################################################################
def parse_header_line(line: str) -> List:
    header_info = json.loads(line[2:])
    if (
        not isinstance(header_info, list)
        or len(header_info) == 0
        or not isinstance(header_info[0], list)
        or len(header_info[0]) != 3
        or header_info[0][0] != "GX taxonomy analysis report"
        or header_info[0][1] < TAXONOMY_ANALYSIS_VER_MIN
        or header_info[0][1] > TAXONOMY_ANALYSIS_VER_MAX
    ):
        return None

    if len(header_info) < 2:
        header_info.append({})

    return header_info


#############################################################################
def action_report(args):

    seqs = []

    # headers
    # ##[["GX taxonomy analysis report",1,1],
    #   {"git-rev":"c2f8891", "run-date":"Mon Jan 31 17:50:42 2022",
    #    "db":{"build-date":"2022-01-28", "seqs":3098, "Gbp":0.59926}}]
    # #seq-id seq-len masked-len cvg-by-all
    #   sep1  tax-name-1  tax-id-1  div-1  cvg-by-div-1  cvg-by-tax-1  score-1
    #   sep2              tax-id-2  div-2  cvg-by-div-2  cvg-by-tax-2  score-2
    #   sep3              tax-id-3  div-3  cvg-by-div-3  cvg-by-tax-3  score-3
    #   sep4              tax-id-4  div-4  cvg-by-div-4  cvg-by-tax-4  score-4
    #   sep5   reserved  result  div      div_pct_cvg
    meta_data_header = ""
    # fields_header = ""
    prev_id = ""
    records = []
    thresholds = ExcludeThresholds(args.good_pct_thr, args.good_len_thr, args.good_abs_len_thr)
    primary_div_name = ""

    with open(args.classify_rpt, "r", encoding="utf8") as fin:
        for line in fin:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith("##"):

                header_info = parse_header_line(line)

                if not header_info:
                    dprint("GX taxonomy analysis report header incorrect")
                    return 1

                header_info[0][0] = "FCS genome report"
                header_info[0][1] = FCS_GENOME_REPORT_VER

                # Read primary div for the file from header
                properties   = header_info[1]
                run_info     = properties.get("run-info", {})
                primary_divs = run_info.get("corrected-primary-divs", run_info.get("primary-divs", []))

                if len(primary_divs) > 0:
                    primary_div_name = primary_divs[0]

                # Fix header for downstream apps, GP-35108
                # remove item with "inferred-primary-divs" from dict run_info even if it does not exist there
                # Revert change in header, GP-35385. Commented line is left as an example
                # run_info.pop("inferred-primary-divs", None)
                # Reconstruct header wirh updated run-info
                header_info[1]["run-info"] = run_info
                meta_data_header = "##" + json.dumps(header_info)

            # elif line[1] == 's':
            #     fields_header = line[1:].rstrip().split("\t")

            if line[0] == "#":
                continue

            r = Record.from_row(line.rstrip().split("\t"))

            if r.orig_id != prev_id:
                # print(f"New id {r.orig_id}", file=sys.stderr)
                if len(records) >= 1:
                    seqs.append(Sequence.calc_sequence(records))

                records = []

            prev_id = r.orig_id
            records.append(r)

        seqs.append(Sequence.calc_sequence(records))

    stats = calc_div_stats(seqs, primary_div_name)
    ungap_len = stats["ungap_len"]
    pa_same_kingdom_contam_demote = stats["pa_same_kingdom_contam_demote"]
    lc_review_divs = stats["lc_review_divs"]

    # do second-pass stuff
    if not primary_div_name:
        primary_div_name = stats["primmary_div_name"]

    if args.debug:
        dprint(f" second_pass check:  {stats}")

    low_cov_mode = False

    if lc_review_divs:
        if args.debug:
            dprint(f"Re-assigning low-coverage for {lc_review_divs} to contam(low-coverage)")

        new_seqs = []

        for seq in seqs:
            new_seqs.append(Sequence.update_low_covs(seq, lc_review_divs))

        seqs = new_seqs
        low_cov_mode = True

    for seq in seqs:
        seq.calc_actions(
            primary_div_name,
            low_cov_mode,
            pa_same_kingdom_contam_demote,
            args.ignore_same_kingdom,
            thresholds
        )

    no_review_without_excludes_rule(seqs)

    if args.production_build_name:
        warning_struct = production_review(seqs, ungap_len, args.production_build_name)

        if warning_struct:
            eprint(json.dumps(warning_struct))

    if meta_data_header != "":
        print(meta_data_header)

    columns = ["seq_id", "start_pos", "end_pos", "seq_len", "action", "div", "agg_cont_cov", "top_tax_name"]

    if args.print_rules:
        columns.append("rules")

    print("#", "\t".join(columns), sep="")

    for seq in seqs:
        for a in seq.actions:
            print(a.tsv(args.print_rules))
    return 0


#############################################################################
def print_rule_codes():
    rules = {
        "G01": "group 'grey' rows with the high-coverage contaminant hits",
        "G02": "group 'grey' rows with the high coverage primary-div hits",
        "G03": "group 'grey' rows with the contaminant, as default",
        "S01": "sequence had only one type of hits: contaminant",
        "S02": "sequence had only one type of hits: primary-div",
        "S03": "sequence had only one type of hits: grey",
        "S04": "sequence is primary-div, putatively extrachromosomal (plastid or mito)",
        "M01": "A possible 'trim' is reduced to exclude for short length and low coverage",
        "M02": "'trim' this range, which is at the sequence edge, of multiple input rows",
        "M03": "'fix' this range, which is in the middle of a sequence, of multiple rows",
        "M04": "'keep' this range, of multiple rows",
        "M05": "'keep' this range, not enough coverage"
    }

    for code, msg in rules.items():
        print(f"  {code}: {msg}")


#############################################################################
def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Classify the outputs of gx-taxify.",
    )

    parser.add_argument(
        "--in",
        required=True,
        dest="classify_rpt",
        help="Output of classify_taxonomy",
    )

# AA: this is not used anywhere
#    parser.add_argument(
#        "--div",
#        dest="div",
#        help="Asserted GX tax-div",
#    )

    parser.add_argument(
        "--rules",
        dest="print_rules",
        action="store_true",
        help="add rules code column",
    )

    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="Extra debug printouts",
    )

    parser.add_argument(
        "--rule-codes",
        dest="print_rule_codes",
        action="store_true",
        help="Print out rule_code list",
    )

    parser.add_argument(
        "--good-pct-thr",
        type=int,
        default=30,
        help="Good sequence percentage threshold",
    )

    parser.add_argument(
        "--good-len-thr",
        type=int,
        default=50000,
        help="Good sequence length threshold",
    )

    parser.add_argument(
        "--good-abs-len-thr",
        type=int,
        default=200,
        help="Good sequence absolute length threshold",
    )

    parser.add_argument(
        "--production-build-name",
        type=str,
        dest="production_build_name",
        help="parameter to set to enable production warnings",
    )

    #added per GP-34669
    parser.add_argument(
        "--ignore-same-kingdom",
        action="store_true",
        dest="ignore_same_kingdom",
        help="parameter to ignore all same kingdom contam, for at least TSA input",
    )

    args = parser.parse_args()

    set_debug(args.debug)

    if args.print_rule_codes is True:
        print_rule_codes()
        return 0

    return action_report(args)


#############################################################################
if __name__ == "__main__":
    sys.exit(main())
