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
# pylint: disable=C0301,C0114,C0103,C0116,C0115,R0914,R0916,R0911
# fmt: off

import sys
assert sys.version_info.major >= 3 and sys.version_info.minor >= 8, f"Python version: {sys.version_info}. Require python 3.8 or newer."

import os
import argparse
import json
import gzip

from dataclasses import dataclass
from collections import defaultdict
from typing      import List, Tuple

# import faulthandler
# import signal
# faulthandler.register(signal.SIGUSR1.value)
#
# To examine current stack-trace without killing the process:
# >>kill -s SIGUSR1 `pidof python`


def as_pct(frac) -> int:  # NB: action-report expects int
    return int(frac * 100 + 0.5)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def is_prok(div) -> bool:
    return div.startswith("prok:") or div.startswith("arch:")


def is_vir_or_synt(div) -> bool:
    # NB: temporary handling for data-error - could be "synt:synthetic" or just "synthetic"
    return div.startswith("virs:") or div.startswith("synt")


def is_euk(div) -> bool:
    return not is_prok(div) and not is_vir_or_synt(div)


vertebrates = [f"anml:{s}" for s in ["mammals", "fishes", "reptiles", "amphibians", "rodents", "birds", "vertebrates", "primates", "marsupials"]]
def tax_kingdom(div): # taxonomic kingdom "prok", "virs", etc; splitting animals into vertabrates and invertebrates
    assert div[4] == ":" or div == "synthetic", div
    return div[:4] if not div.startswith("anml:") else "vrbt" if div in vertebrates else "ivrt"


def is_same_kingdom(div1: str, div2: str) -> bool:
    null_divs = ["NULL", "synthetic"]
    if div1 in null_divs or div2 in null_divs:
        return False
    assert div1[4] == ":" and div2[4] == ":", (div1, div2)
    return tax_kingdom(div1) == tax_kingdom(div2) and (not tax_kingdom(div1) == "anml" or (div1 in vertebrates) == (div2 in vertebrates))

assert     is_same_kingdom("anml:mammals", "anml:fishes")
assert not is_same_kingdom("anml:mammals", "anml:nematodes")  # vertebrate / not-vertebrate
assert not is_same_kingdom("anml:mammals", "prok:bacteria")


def as_readable(num) -> str:
    s = "PTGMK "
    while len(s) > 1 and abs(num) > 1000:
        (num, s) = (num / 1000, s[:-1])
    return f"{num} bp" if s[-1] == " " else f"{int(num * 100 + 0.5) / 100} {s[-1]}bp"

assert as_readable(123456) == "123.46 Kbp"


#############################################################################
@dataclass
class Taxon:
    # fmt: off
    div        : str = ""  # tax-div, e.g. anml:primates. Prefix is the tax-kingdom
    tax_id     : int = 0   # e.g. 30033
    cvg_by_div : int = 0   # absolute coverage (length) by all alignments in the same div
    cvg_by_tax : int = 0   # absolute coverage (length) by alignments of this taxon
    score      : int = 0   # aggregate score of the alignments of this taxon.
    # fmt: on

    @staticmethod
    def from_row(row: List[str], i: int):
        assert i in [6, 12, 18, 24]  # starting position for taxon-subrecord
        return Taxon(row[i + 1], *(int(row[i + k]) for k in (0, 2, 3, 4)))


#############################################################################
# Corresponds to a single row in gx-taxify output
@dataclass
class Record:
    # fmt: off
    seq_id          : str  # ~ indicates subseq split on on N-runs; ~~ subseq split on chimeric switchpoint.
    len             : int  # sequence length
    transposon_len  : int  # transposons coverage length
    repeat_len      : int  # low-complexity + conserved regions length
    cvg_by_all      : int  # coverage length by all alignments
    taxa            : List[Taxon] # len up-to 4, may be empty
    # fmt: on

    #########################################################################
    @staticmethod
    def from_row(row: List[str]):
        """Create Record from fields of taxonomy report file (output of gx-taxify)"""

        assert len(row) >= 30 and row[4] == "|" and row[29] == "|"

        repeats_lens = [int(s) for s in row[2].split(",")]  # [transposons, low-complexity, conserved, n-runs]
        assert len(repeats_lens) == 4

        return Record(
            row[0],                             # seq-id
            int(row[1]) - repeats_lens[3],      # seq-len excluding Ns
            repeats_lens[0],                    # transposons
            repeats_lens[1] + repeats_lens[2],  # low_complexity + conserved
            int(row[3]),                        # cvg-by-all
            [Taxon.from_row(row, i) for i in [6, 12, 18, 24] if row[i + 1]]
        )

    #########################################################################
    @staticmethod
    def get_rows(taxonomy_rpt: str):  # -> sequence of (metaline or list-of-fields)
        with (
            gzip.open(taxonomy_rpt, "rt", encoding="utf8") if taxonomy_rpt.endswith(".gz")
            else open(taxonomy_rpt, "rt", encoding="utf8")
        ) as fin:
            for line in fin:
                if line:
                    yield line if line.startswith("#") else line.rstrip().split("\t")[:30]


    #########################################################################
    def get_top_taxon_frac_coverage(self) -> float:
        assert self.taxa
        t0 = self.taxa[0]
        eff_len = max(self.len - self.repeat_len - self.transposon_len, self.cvg_by_all)
        return t0.cvg_by_div / eff_len

    #########################################################################
    def is_top_taxon_high_evidence(self) -> bool:
        # Distributions of score^2 vs. length are different for same-div vs spurious alignments, and they are
        # linearly separable.
        #
        # Collect (len, score) pairs:
        # >> cat GCF_002007445.1.taxonomy.rpt | tail -n+3 | awkt '1{ print $10,$11; print $16,$17; print $22,$23; print $28,$29 }' \
        #                                                 | awkt '($1 > 50){ print $1, $2*$2 - $1*15 }' > tmp.pts;
        # >> gnuplot -e 'plot "tmp.pts" with dots; pause -1'
        #
        # The upper-arm (positive y-values) comes from same-div alignments; the lower arm is from spurious hits.
        t0 = self.taxa[0]
        is_high_scoring = t0.cvg_by_tax > 200 and t0.score * t0.score > 15 * t0.cvg_by_tax
        return is_high_scoring and self.get_top_taxon_frac_coverage() > 0.5


#############################################################################
def get_metadata(taxonomy_rpt: str) -> list:
    ret = next((json.loads(row[2:]) for row in Record.get_rows(taxonomy_rpt) if isinstance(row, str) and row.startswith("##")), None)

    assert isinstance(ret, list) and len(ret) > 1 and isinstance(ret[0], list) and len(ret[0]) == 3  # [format-name, major-ver, minor-ver]
    assert ret[0][0] in ("GX taxonomy pre-analysis report", "GX taxonomy analysis report") and ret[0][1] <= 3

    return ret

#############################################################################
def classify_record(
    r             : Record,
    primary_divs  : List[str],
    contam_divs   : List[str],
    min_cvg_frac  : float,
    is_outlier_org: bool = False,
) -> Tuple[str, str]:  # (div, class-label), e.g. ("anml:fishes", "primary-div")

    # NB: each return below has different class-label, so that a result can be tracked to corresponding rule.
    # Possible class-labels (see GP-34548):
    # Primaries     - primary-div, repeat, transposon
    # Contaminants  - contaminant(div), contaminant(virus), contaminant(cross-kingdom),
    #                 contaminant(cross-div), contaminant(human), same-kingdom-chimeric
    # Inconclusives - low-coverage, inconclusive, inconclusive(metagenome)

    is_metagenome = primary_divs[0] == "unkn:metagenomes"
    is_low_cvg    = not r.taxa or r.taxa[0].cvg_by_div < min_cvg_frac * r.len
    is_high_cvg   = r.taxa and r.taxa[0].cvg_by_div > 0.8 * r.len
    is_transposon = r.transposon_len > r.taxa[0].cvg_by_div if r.taxa else 0.3 * r.len
    is_repeat     = r.repeat_len + r.transposon_len > min(r.cvg_by_all, 0.66 * r.len) # GP-31287

    if is_low_cvg or is_transposon:
        label = "repeat" if is_repeat else "low-coverage" if is_low_cvg else "transposon"
        return (r.taxa[0].div if r.taxa else "none", label)

    # div's preference-facter when selecting candidate divs and taxa below
    def get_factor(div):
        return 0.5 if div in primary_divs or is_vir_or_synt(div) else 0.8

    # Filter taxa to entries that are close to t0 by product of score and cvg_by_div, more leeway for primary-div taxa.
    # unit-test cases 5,6
    t0 = r.taxa[0]
    taxa = [t for t in r.taxa if t.score * t.cvg_by_div > t0.score * t0.cvg_by_div * get_factor(t.div)]
    divs = [t.div for t in taxa]
    kdms = list(set((tax_kingdom(t.div) for t in taxa)))
    is_same_kdm = tax_kingdom(primary_divs[0]) in kdms

    # GP-31427
    # Special handling for human contamination even it's in primary-divs, e.g. bat GCA_014824575.1, panda GCF_002007445.1
    # Need to be conservative, because some species may be half-way between rodents and primates, e.g. mole rat GCA_000247695.1
    #
    # NB: this rule goes before "primary-div" calls.
    if (t0.div == "anml:primates"
        and primary_divs[0] != "anml:primates"
        and any((t.tax_id == 9606 for t in r.taxa[:2]))
        and r.get_top_taxon_frac_coverage() > 0.7
        and r.is_top_taxon_high_evidence()
        and "~~" not in r.seq_id  # not chimeric
        and not any((t.score > t0.score * 0.8 and t.div != t0.div for t in taxa))
    ):
        return (t0.div, "contaminant(human)")

    for div in divs:  # primary-div calls.
        if (   div in primary_divs
            or is_metagenome and (is_prok(div) or is_vir_or_synt(div))
            or is_same_kdm and is_outlier_org  # Treat same-kingdom hits as primary-divs if outlier-org. GP-34615
        ):
            return (div, "primary-div")

    for div in divs:  # virus calls. GP-34699, GP-34839
        if div.startswith("virs:"):
            is_prokv_in_prok = is_prok(primary_divs[0]) and div in ["virs:viruses", "virs:prokaryotic viruses"]
            is_eukv_in_euk   =  is_euk(primary_divs[0]) and div in ["virs:viruses", "virs:eukaryotic viruses"]
            is_integrant     = "~" in r.seq_id  # either contig-level ("~") or chimeric ("~~")
            treat_as_primary = is_prokv_in_prok or (is_eukv_in_euk and is_integrant)

            return (div, "primary-div(virus)") if treat_as_primary else (div, "contaminant(virus)")

    # same-kingdom chimeric integrant are likely FPs. GCA_918807975.1 (tunicate) has lots of these.
    if ("~~" in r.seq_id and is_same_kdm and r.len < 10000):  # pylint: disable=R1705
        return (t0.div, "same-kingdom-chimeric")

    # NB: GCA_900118605.1 (chlamydia) has cases of nematode-or-rodent calls (len(kdms)==2) that are both in contam_divs
    elif (((t0.div in contam_divs and len(kdms) == 1) or all((d in contam_divs for d in divs)))
        and (t0.score > 50 or (is_high_cvg and not is_same_kdm))
    ):
        return (t0.div, "contaminant(div)")

    elif is_metagenome:
        # plants, fungi, protists, invertebrates, low-scoring contam-vertebrates.
        # Should "large" invertebrates go on the contam_divs too?
        return (t0.div, "inconclusive(metagenome)")

    elif (  len(kdms) == 1
        and not is_same_kdm
        and t0.score > 100
        and (r.is_top_taxon_high_evidence())
    ):
        return (t0.div, "contaminant(cross-kingdom)")

    elif (  len(divs) == 1
        and t0.div not in primary_divs
        and t0.score > 100
        and (r.is_top_taxon_high_evidence() or sum((1 for t in taxa if t.div == t0.div)) >= 2)
        and (not is_repeat or t0.div == "anml:primates")
    ):
        return (t0.div, "contaminant(cross-div)")
    else:
        return (t0.div, "repeat" if is_repeat else "inconclusive")


#############################################################################
def select_divs(taxonomy_rpt: str, primary_div: str) -> List[str]:
    """ return list of divs that are "well-represented" in taxonomy.rpt """

    d = defaultdict(lambda: [0, 0, 0, 0, 0])  # div -> vec-of-aggregates; see directly below for what these contain.

    for r in (Record.from_row(row) for row in Record.get_rows(taxonomy_rpt) if isinstance(row, list)):
        cvg_by_primary_div      = max((t.cvg_by_div for t in r.taxa if t.div == primary_div), default=0)
        for t in r.taxa:
            # For better SNR ignore low-scoring noise.
            # However, for highly fragmented assemblies such as GCA_001891065.1
            # even the threshold value of 150 is too high, so alternatively take if high-coverage.
            if t.score > r.taxa[0].score * 0.8 and (t.score > 150 or t.cvg_by_tax > r.len * 0.8):
                d[t.div][0] += t.cvg_by_div
                d[t.div][1] += min(t.cvg_by_div, r.repeat_len)  # cvg-by-div possibly repeat-specific
                d[t.div][2] += min(t.cvg_by_div, r.transposon_len)  # cvg-by-div possibly transposon-specific
                d[t.div][3] += min(t.cvg_by_div, cvg_by_primary_div)  # cvg-by-div possibly primary-div-specific
                d[t.div][4] += max(0, t.cvg_by_div - r.repeat_len - r.transposon_len) # cvg-by-div not repeat-or-transposon-specific

    items = sorted(d.items(), key=lambda item: item[1][0], reverse=True)  # by decreasing sum cvg_by_div

    # divs that have excess of aggregate cvg_by_div vs. repeats-or-transposon of >10kb (last field in above vectors),
    # and the divs are not exceedingly repeat-specific or transposon-specific or primary-div-specific,
    # e.g. algae in fish GCA_001640805.2, or marsupials in armadillo GCA_000208655.2
    def select(v):
        return max(v[1], v[2], v[3]) < v[0] * 0.75 and v[4] > 10000

    if os.getenv("GX_CLASSIFY_TAXONOMY_VERBOSE") and items:
        eprint("\nTop represented putative divs:")
        eprint("#\t", "coverage", "repeats_pct", "transposon_pct", "primary_div_pct", "div")
        for div, v in items[:10]:
            eprint("\t", as_readable(v[0]), *(as_pct(v[i] / v[0]) for i in (1, 2, 3)), "T" if select(v) else "F", div, sep="\t")
        eprint("")

    return [div for (div, v) in items if div == items[0][0] or select(v)]


#############################################################################
def adjust_divs(asserted_div, inferred_primary_divs): # -> (primary_divs, contam_divs)
    contam_divs = []
    primary_divs = inferred_primary_divs.copy()

    # GP-33376
    # prok and fung divs having >100 assemblies in the db.
    # cat /dev/shm/gxdb/all.assemblies.tsv | awkt '$6' | cut -f9 | tabulate | grep -P 'prok:|fung:' | awkt '($2 > 100) {print "\""$1"\", " }'
    well_represented_divs = [
        "prok:firmicutes",
        "prok:high GC Gram+",
        "prok:g-proteobacteria",
        "prok:a-proteobacteria",
        "prok:CFB group bacteria",
        "fung:ascomycetes",
        "prok:b-proteobacteria",
        "fung:basidiomycetes",
        "fung:budding yeasts",
        "prok:cyanobacteria",
        "prok:d-proteobacteria",
        "prok:proteobacteria",
        "prok:actinobacteria",
        "prok:mycoplasmas",
        "prok:spirochetes",
    #   "prok:bacteria",  # excluding generic divs - see GP-33376
    #   "fung:fungi",
    ]

    # Replace-or-append primary-divs with asserted-div as appropriate.
    if asserted_div == "virs:viruses":  # GP-33387, GP-34316
        primary_divs = ["virs:viruses", "virs:eukaryotic viruses", "virs:prokaryotic viruses"]

    elif asserted_div in ["virs:eukaryotic viruses", "virs:prokaryotic viruses"]:
        primary_divs.append("virs:viruses")

    elif asserted_div == "unkn:metagenomes": # GP-33795
        primary_divs = ["unkn:metagenomes"]
        contam_divs = vertebrates  # should be chordates, but close enough. What about "large" invertebrates?

    elif asserted_div is None or asserted_div in primary_divs:
        pass  # Normal case

    elif is_same_kingdom(asserted_div, primary_divs[0]) and asserted_div not in well_represented_divs:
        primary_divs = [asserted_div] + primary_divs
        eprint(f"\n    * * * Adding asserted tax-div '{asserted_div}' to the set of primary divs. * * *\n")

    else :  # GP-33376: if different-kingdom or asserted-div is high-confidence but missing from primary-divs,
            # then trust the asserted-div, and treat the inferred-divs as contamination (similar to metagenomes case).
        primary_divs = [asserted_div]
        contam_divs = inferred_primary_divs.copy()

        eprint("\n\n--------------------------------------------------------------------------------------------------")
        eprint(f"Warning: Asserted tax-div '{asserted_div}' is " + (
                "from a different tax-kingdom than the inferred-primary-divs." if not is_same_kingdom(asserted_div, primary_divs[0])
           else "well-represented in db, but absent from inferred-primary-divs."
        ))
        eprint("This means that either asserted tax-div is incorrect, or the input is predominantly contamination.")
        eprint("Will trust the asserted div and treat inferred-primary-divs as contaminants.")
        eprint("--------------------------------------------------------------------------------------------------\n")

    return primary_divs, contam_divs


#############################################################################
def classify_taxonomy(args):

    metadata                  = get_metadata(args.taxonomy_rpt)
    run_info                  = metadata[1]["run-info"]
    agg_cvg_frac              = run_info["agg-cvg"]
    inferred_primary_divs     = run_info["primary-divs" if "primary-divs" in run_info else "inferred-primary-divs"]
    asserted_div              = run_info.get("asserted-div", None)
    primary_divs, contam_divs = adjust_divs(asserted_div, inferred_primary_divs)
    selected_divs             = select_divs(args.taxonomy_rpt, primary_divs[0]) if asserted_div != "unkn:metagenomes" else []
    top_div                   = selected_divs[0] if selected_divs else primary_divs[0]

    is_outlier_org = ( # GP-34615 - egregious contamination or "weird-bastie" - will treat same-kingdom divs as primary-divs.
        asserted_div
        and asserted_div not in inferred_primary_divs
        and not is_same_kingdom(asserted_div, inferred_primary_divs[0])
    )

    #########################################################################
    if top_div not in primary_divs and top_div not in contam_divs and is_same_kingdom(top_div, primary_divs[0]):
        # E.g. Alitta segmented worm - asserted-div is worms, primary-div by top-cvg from taxify is insects; top-div is molluscs
        primary_divs.append(top_div)
        eprint(f"\n   * * * Adding dominant div '{top_div}' to the set of primary divs. * * * \n")

    contam_divs += [div for div in selected_divs if div not in primary_divs and div not in contam_divs]
    contam_divs += [s.strip() for s in os.environ.get("GX_EXTRA_CONTAM_DIVS", "").split(",") if s]
    assert all(div[4] == ":" or div in ["synthetic"] for div in contam_divs), contam_divs

    #########################################################################
    eprint("Asserted div               :", asserted_div)
    eprint("Inferred primary-divs      :", inferred_primary_divs)
    eprint("Corrected primary-divs     :", primary_divs)
    eprint("Putative contaminant divs  :", contam_divs)
    eprint(f"Aggregate coverage         : { int(agg_cvg_frac * 100 + 0.5) }%")

    # For high-coverage genomes a no-coverage sequence is indicative of novel contaminant,
    # whereas for low-coverage genomes it's indicative of primary-div - hence require higher minimum coverage.
    min_cvg_frac = max(0.2, 0.6 * (1 - agg_cvg_frac))
    eprint(f"Minimum contam. coverage   : { as_pct(min_cvg_frac) }%")

    out_filepath = None if not args.out_dir else args.out_dir + "/" + os.path.basename(args.taxonomy_rpt)
    fout = (
        sys.stdout if not args.out_dir
        else gzip.open(out_filepath, "wt", encoding="utf8") if out_filepath.endswith(".gz")
        else      open(out_filepath, "wt", encoding="utf8")  # pylint: disable=R1732
    )

    for row in Record.get_rows(args.taxonomy_rpt):
        if isinstance(row, str) and row.startswith("##"):
            assert any(s in row for s in ("GX taxonomy pre-analysis report", "GX taxonomy analysis report"))
            assert "run-info" in metadata[1]
            metadata[0][0] = "GX taxonomy analysis report"
            metadata[1]["run-info"]["corrected-primary-divs"] = primary_divs  # GP-34560
            print("##" + json.dumps(metadata), file=fout)
        elif isinstance(row, str) and row.startswith("#"):
            print(row.rstrip(), "reserved", "result", "div", "div_pct_cvg", sep="\t", file=fout)
        else:
            r = Record.from_row(row)
            (div, label) = classify_record(r, primary_divs, contam_divs, min_cvg_frac, is_outlier_org)
            div_cvg = next((t.cvg_by_div for t in r.taxa if t.div == div), 0)
            # NB: r.len can be 0, e.g. all-Ns, GP-34773
            print(*row, "n/a", label, div, as_pct(div_cvg / r.len if r.len else 0), sep="\t", file=fout)

    fout.close()

#############################################################################
def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Classify the outputs of gx-taxify.",
    )

    parser.add_argument("--in", dest="taxonomy_rpt", type=str, required=True, help="Output of gx taxify")
    parser.add_argument("--out-dir", type=str, required=False, help="If specified, put the output in the file with the same basename in --out-dir; otherwise output to stdout.")
    classify_taxonomy(parser.parse_args())
    return 0


#############################################################################
if __name__ == "__main__":
    sys.exit(main())
