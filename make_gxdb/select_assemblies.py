#!/usr/bin/env python3
# pylint: disable=C0301,C0114,C0103,C0116,R0913,W0621,R0914
# fmt: off

# Prepare assemblies.tsv and related files that are used
# by make_gxdb.sh to create GX database.
#
# Here we select a minimal, and diverse set of assemblies
# from RefSeq and GenBank.

import argparse
import sys
import os
import gzip
import subprocess

from collections import defaultdict, namedtuple, Counter
from xml.etree import cElementTree as ET
from pathlib import Path
from typing import Dict, List, Set, NewType, Any


# fmt: off
TaxId    = NewType('TaxId',    int)
GcAccVer = NewType('GcAccVer', str)  # GenColl GCA or GCF accession.version
BlastDiv = NewType('BlastDiv', str)  # e.g. "primates"
GxDiv    = NewType('GxDiv',    str)  # e.g. "anml:primates"
# fmt: on

#############################################################################


def eprint(*args):
    print(*args, file=sys.stderr, flush=True)


scriptdir = os.path.realpath(
    os.path.dirname(os.readlink(sys.argv[0]) if Path(sys.argv[0]).is_symlink() else sys.argv[0])
)

#############################################################################

parser = argparse.ArgumentParser(
    description="Select assemblies in scope for gxdb.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

# GP-33261
# TSV of tax-id-1, tax-id-2, similarity-score (higher is better).
parser.add_argument(
    "--tax-distances",
    default=os.getenv("FCS_ROOT", ".") + "/gx/data/taxa_proximities/taxa.proximities.tsv.gz",
    help="taxonomic proximities file.",
)

parser.add_argument(
    "--gencoll-sqsh-ms-args",
    default=os.getenv("GENCOLL_DB_AUTH", "-S server-name -D db-name -U login -P password"),
    help="sqsh-ms connection args for gencoll-db.",
)

parser.add_argument(
    "--gc-genomes-root-dir",
    default=os.getenv("GC_GENOMES_ROOT_DIR"),
    help="Local path to https://ftp.ncbi.nlm.nih.gov/genomes/",
)

parser.add_argument(
    "--out-dir",
    default="./assemblies_out",
    help="Output directory.",
)

args = parser.parse_args()

#############################################################################

if not Path(args.out_dir).exists():
    os.mkdir(args.out_dir)

os.chdir(args.out_dir)
eprint(f"Output directory: {args.out_dir}")


#############################################################################
def prev_accver(accver):
    if "." not in accver:
        return accver
    acc, ver = accver.split(".")
    ver = int(ver) - 1
    return f"{acc}.{ver}"


assert prev_accver("GCF_000001405.40") == "GCF_000001405.39"


#############################################################################
def run(cmd: str):
    cmd = " ".join(cmd.split())  # collapse newlines and multiple spaces
    eprint("\nExecuting: ", cmd, "\n...")
    subprocess.run(f"set -eo pipefail; {cmd}", shell=True, check=True)


#############################################################################
def generate_if_not_exists(cmds, out_filename: str):
    """
    Run commands to produce the file unless it already exists.

    (Most of our commands to fetch the data take a while,
    so we'll dump their output into files, and read it from the
    files on afterwards, and on re-executiion.)
    """
    if Path(out_filename).exists():
        eprint(f"Reusing existing {out_filename}.")
    else:
        for cmd in (cmds if isinstance(cmds, list) else [cmds]):
            run(cmd)

    assert Path(out_filename).exists()


#############################################################################
def load_dict(filename: str, fn) -> Dict[Any, Any]:
    """
    Populate a dictionary from a TSV file.
    fn is a function taking cols:List[str], returning key-value pair or None.
    """
    d = {}
    eprint(f"Loading {filename}...")
    with (
        gzip.open(filename, "rt", encoding="utf8")
        if filename.endswith(".gz")
        else open(filename, "rt", encoding="utf8")
    ) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            kv = fn(line.rstrip("\n").split("\t"))
            if kv is None:
                continue
            val = d.get(kv[0])
            if val is None:
                d[kv[0]] = kv[1]
            elif val != kv[1]:
                eprint(f"Warning: Key-conflict in {filename}: {kv}; existing-value: {val}")
    return d


#############################################################################
def generate_dict(cmds, filename: str, fn) -> Dict[Any, Any]:
    generate_if_not_exists(cmds, filename)
    return load_dict(filename, fn)


#############################################################################
tax2blast_div: Dict[TaxId, BlastDiv] = generate_dict(
    "taxpropdump -t $'\\t' blast | sed 's/\"//g' | gzip -c > taxpropdump.tsv.gz",
    "taxpropdump.tsv.gz",
    lambda cols: (TaxId(int(cols[0])), BlastDiv(cols[1])),
)


blast_div2gx: Dict[BlastDiv, GxDiv] = load_dict(
    f"{scriptdir}/blast_names_mapping.tsv",
    lambda cols: (BlastDiv(cols[0]), GxDiv(cols[1])),
)

# Category is <RefSeq_category> from esummary xml,
# or fifth column in assembly_summary_refseq.txt.
#
# Will override refseq-category by GCF from assembly_summary_refseq.txt,
# to prefer GCF reference assemblies where there's a newer reference genbank
# without corresponding refseq, e.g. prefer GCA_000001405.28/GCF_000001405.39
# over latest GCA_000001405.29.
#
# <RefSeq_category/> element is 'na' for GCF_000001405.39 in the assembly docsum,
# whereas "reference genome" in assembly_summary_refseq.txt.
gc2category: Dict[GcAccVer, str] = load_dict(
    args.gc_genomes_root_dir + "/ASSEMBLY_REPORTS/assembly_summary_refseq.txt",
    lambda cols: (GcAccVer(cols[0]), cols[4]),
)

#############################################################################
# Will use it to determine "taxon-popularity" when sorting assemblies by preference.
#
# Alternatively, can collect from docsums and populate Dict[TaxId, Set[str]]
# xtract ... -block GB_BioProjects -element Bioproj/BioprojectAccn
tax2num_bioprojects: Dict[TaxId, int] = generate_dict(
    f""" cat {args.gc_genomes_root_dir}/ASSEMBLY_REPORTS/assembly_summary_genbank.txt |
         tail -n+3 |
         datamash -s -g6 countunique 2 > bioproject_counts.tsv
    """,
    "bioproject_counts.tsv",
    lambda cols: (TaxId(int(cols[0])), int(cols[1])),
)

#############################################################################
# Will distinguish prok and euk viruses.
# There are two ways to go about it: infer via translation-table,
# or trust prokaryotic-virus[orgn].
# These methods overlap by >98%. We'll take the union of both.
prok_vir_taxa: Set[TaxId] = set(
    generate_dict(
        [
            """ esearch -db taxonomy -query 'viruses[orgn] translation-table-11[prop]' |
                efetch -format uid > prok_viruses.tsv.1
            """,
            """ esearch -db taxonomy -query 'viruses[orgn] translation-table-11[prop]' |
                elink -related -name taxonomy_taxonomy_child |
                efetch -format uid > prok_viruses.tsv.2
            """,
            """ esearch -db taxonomy -query 'prokaryotic-virus[orgn]' |
                efetch -format uid > prok_viruses.tsv.3
            """,
            "cat prok_viruses.tsv.1 prok_viruses.tsv.2 prok_viruses.tsv.3 | sort -u > prok_viruses.taxids",
            "rm prok_viruses.tsv.1 prok_viruses.tsv.2 prok_viruses.tsv.3",
        ],
        "prok_viruses.taxids",
        lambda cols: (TaxId(int(cols[0])), True),
    ).keys()
)


# NB 'viruses[orgn] translation-table-1[prop]' yields a superset of translation-table-11 used to get prok-taxa - can't use it.
euk_vir_taxa: Set[TaxId] = set(
    generate_dict(
        [
            """ esearch -db taxonomy -query 'viruses[Organism] AND (vhost-algae[Filter] OR vhost-fungi[Filter] OR vhost-invertebrates[Filter] OR vhost-plants[Filter] OR vhost-protozoa[Filter] OR vhost-vertebrates[Filter] OR vhost-eukaryotes[Filter]) NOT (archaea[Organism] OR bacteria[Organism])' |
                efetch -format uid > euk_viruses.taxids
            """
        ],
        "euk_viruses.taxids",
        lambda cols: (TaxId(int(cols[0])), True),
    ).keys()
)


#############################################################################
generate_if_not_exists(
    [
        f"sqsh-ms {args.gencoll_sqsh_ms_args} -i {scriptdir}/V_ContamGX.sql",
        "touch V_ContamGX.done",
    ],
    "V_ContamGX.done",
)

# The above comand produces these three files:
assert Path("gc_screened_assemblies.tsv").exists()
assert Path("gc_assemblies_denylist.tsv").exists()
assert Path("gc_sequences_denylist.tsv").exists()

screened_asms: Set[GcAccVer] = set(
    load_dict("gc_screened_assemblies.tsv", lambda cols: (GcAccVer(cols[0]), True)).keys(),
)


#############################################################################
def get_asm_control(tax_id: int, gca: GcAccVer, gcf: GcAccVer) -> (bool, bool):
    """return (is-on-acceptlist, is-on-denylist)"""

    # load .asm_dict and .tax_dict on the first invocation
    if not hasattr(get_asm_control, "asm_dict"):
        get_asm_control.asm_dict = load_dict(
            "gc_assemblies_denylist.tsv",
            lambda cols: (GcAccVer(cols[2]), cols[0]),
        )

        # Augment/override from curated assemblies_control.tsv:
        get_asm_control.asm_dict.update(
            load_dict(
                f"{scriptdir}/assemblies_control.tsv",
                lambda cols: (GcAccVer(cols[2]), cols[0]),
            )
        )

        # To accept/reject by tax-id rather than assembly-id:
        get_asm_control.tax_dict = load_dict(
            f"{scriptdir}/assemblies_control.tsv",
            lambda cols: None if cols[2] != "*" else (TaxId(int(cols[1])), cols[0]),
        )

    assert tax_id > 0
    assert gca.startswith("GCA_") or not gca
    assert gcf.startswith("GCF_") or not gcf

    val = get_asm_control.asm_dict.get(gcf, get_asm_control.asm_dict.get(gca))
    assert val in (None, "+", "-")  # '+' means on-acceptlist; '-' means on-denylist.

    is_accept = val == "+"
    is_reject = val == "-" or (get_asm_control.tax_dict.get(tax_id) == "-" and not is_accept)
    return (is_accept, is_reject)


#############################################################################
def fetch_assemblies_esummaries(out_filename):

    if Path(out_filename).is_file():
        eprint(f"Reusing existing {out_filename}.")
        return

    # GP-32020
    mag_exclusions = [
        "abnormal gene to sequence ratio",
        "chimeric",
        "contaminated",
        "derived from environmental source",
    #   "derived from single cell",
        "fragmented assembly",
        "from large multi isolate project",
        "genome length too large",
        "genome length too small",
    #   "genus undefined",
        "hybrid",
        "low gene count",
        "low quality sequence",
        "many frameshifted proteins",
        "metagenome",
        "misassembled",
        "missing ribosomal protein genes",
        "missing rrna genes",
    #   "missing strain identifier",
        "missing trna genes",
        "mixed culture",
        "not used as type",
        "partial",
        "refseq annotation failed",
    #   "sequence duplications",
        "unverified source organism",
    ]

    mag_exclusions_str = " OR ".join((f'"{ex}"[Excluded from RefSeq]' for ex in mag_exclusions))
    prok_mags = f'(prokaryotes[orgn] AND "derived from metagenome"[Excluded from RefSeq] NOT ({mag_exclusions_str}))'

    assembly_query = (
        "      (latest_genbank[filter] OR latest_refseq[filter])"
        f" AND (representative[prop] OR eukaryotes[orgn] OR viruses[orgn] OR {prok_mags})"
        "  NOT (replaced_refseq[prop] OR suppressed_refseq[prop] OR anomalous[filter])"
    )

    run(f"esearch -db assembly -query '{assembly_query}' > {out_filename}.esearch")
    num_assemblies_from_esearch = int(ET.parse(f"{out_filename}.esearch").getroot().find("Count").text)
    eprint("Total assemblies form esearch:", num_assemblies_from_esearch)

    # NB: We mangle the Stats-XML with awk below so that Meta/Stat element contains total_length
    # ExclFromRefSeq/string is variable number of columns, so it's placed last.
    #
    # NB: can't use empty-string as NULL placeholder, as it will swallow columns.
    run(
        f"""
    cat {out_filename}.esearch |
        esummary |
        awk '$0 !~ /Stats>/ && ($0 !~ /<Stat category/ || $0 ~ /<Stat category="total_length"/)' |
        xtract -pattern DocumentSummary           \
                    -element Taxid,Meta/Stat,ContigN50,SpeciesTaxid,SpeciesName,Organism \
          -def NULL -element Synonym/Genbank      \
          -def NULL -element Synonym/RefSeq       \
          -def NULL -element AssemblyName         \
          -def NULL -element Id                   \
          -def NULL -element RefSeq_category      \
          -def NULL -element FtpPath_Assembly_rpt \
          -def NULL -element ExclFromRefSeq/string |
        pv --line-mode -bat > {out_filename}.tmp
    """
    )
    os.rename(f"{out_filename}.tmp", out_filename)
    os.remove(f"{out_filename}.esearch")

    # Sometimes esummary has internal failures and returns partial output without
    # erroring-out, so we must manually verify the number of returned entries.
    with open(out_filename, encoding="ascii") as f:
        assert num_assemblies_from_esearch == sum(1 for _ in f)


#############################################################################
Assembly = namedtuple(
    "Assembly",
    "tax_id len n50 species orgname gca gcf blast_div gx_div "
    "asm_name asm_id rs_category excl_from_refseq weight num_bioprojects sp_tax_id nbr_tax_id status ftp_path",
)

genus_div_override = dict({
    "Tupaia"      : "anml:rodents",
    "Galeopterus" : "anml:primates",
    "Alligator"   : "anml:reptiles",
    "Gavialis"    : "anml:reptiles",
    "Manis"       : "anml:mammals",
})


# Line is one record from fetch_assemblies_esummaries() report above
# (output of the xtract command).
def line2asm(line: str) -> Assembly:

    #########################################################################
    # Fields from xtract DocumentSummary command above.

    # Convert "NULL"s to empty-strings, so we can rely on conversion-to-bool.
    cols = ["" if s == "NULL" else s for s in line.split("\t")]

    (tax_id, asm_len, n50, sp_tax_id) = (int(x) for x in cols[0:4])
    (species, orgname, gca, gcf, asm_name, asm_id, rs_category, ftp_path) = cols[4:12]
    asm_id = int(asm_id)  # type: ignore
    excl_from_refseq = [x for x in cols[12:] if x]  # may be more than one.

    assert gca.startswith("GCA_") or not gca
    assert gcf.startswith("GCF_") or not gcf

    assert ftp_path.startswith("ftp://") or not ftp_path

    rs_category = gc2category.get(gcf, rs_category)
    assert rs_category in ("na", "representative genome", "reference genome")

    #########################################################################

    # orgname is e.g. "Chitinophaga pinensis DSM 2588 (CFB group bacteria)" or "Homo sapiens (human)"
    assert "(" in orgname and orgname[-1] == ")"
    name_in_paren = orgname[orgname.find("(") + 1 : -1]  # extrac blast-div or common-name from parens.

    # Compute blast-div.
    blast_div = tax2blast_div.get(tax_id) or (name_in_paren if name_in_paren in blast_div2gx else "unknown")
    if blast_div == "unknown":
        eprint("Could not look up blast-div for", tax_id)

    # Compute gx-div and kingdom.
    gx_div = blast_div2gx[blast_div]
    assert gx_div and gx_div[4] == ":"
    gx_kingdom = gx_div[0:4]
    if gx_kingdom == "virs":
        gx_div = (
                 "virs:prokaryotic viruses" if " phage " in species or species.endswith(" phage") # GP-34479
            else "virs:prokaryotic viruses" if tax_id in prok_vir_taxa and tax_id not in euk_vir_taxa
            else "virs:eukaryotic viruses"  if tax_id in euk_vir_taxa and tax_id not in prok_vir_taxa
            else "virs:viruses"
        )

    if gx_kingdom not in ("virs", "prok"):
        gx_div = genus_div_override.get(species.split(" ")[0], gx_div)
        assert gx_kingdom == gx_div[0:4]


    # Replace orgname with the name from paretheses, unless it's div.
    if name_in_paren != blast_div:
        orgname = name_in_paren

    (is_accept, is_reject) = get_asm_control(tax_id, gca, gcf)

    if sp_tax_id == 10847:  # GP-34660 - excluding phiX viruses; will be included as synthetic via gcontam1.
        is_reject = True

    # fmt: off

    # Will use it to determine the assemblies preference order.
    def get_weight():
        is_repr = rs_category.startswith("representative")
        is_ref = rs_category.startswith("reference")

        # taxa in scope for GenBank supplementation:
        gb_ok = (
            gx_div in ("anml:rotifers", "anml:nematodes", "plnt:green algae")
            or gx_kingdom in ("fung", "prst")

            # MAGs. for now excluding MAGs that have extra exclusions, e.g. "genus undefined",
            # or gx_kingdom == "prok" and "derived from metagenome" in excl_from_refseq # MAGs
            or gx_kingdom == "prok" and excl_from_refseq == ["derived from metagenome"]

            or gx_kingdom == "virs" and not excl_from_refseq  # genbank viruses
        )

        return (
            -1      if is_reject
            else 10 if is_accept
            else 7  if gcf and is_ref
            else 6  if gcf and is_repr
            else 5  if gcf
            else 4  if gb_ok and is_ref
            else 3  if gb_ok and is_repr
            else 2  if gb_ok
            else 1  if is_repr or is_ref
            else 0
        )

    weight = get_weight()

    # GP-32813, GP-32813
    is_prescreened = (
        blast_div == "viruses"  # these are not in scope for prescreening
        or weight >= 7  # reference or on-acceptlist
        or (gca and gca in screened_asms)
        or (gcf and gcf in screened_asms)
        or (prev_accver(gcf) in screened_asms and gx_div[0:4] not in ("prok", "arch"))
    )

    #eprint("is_prescreened = True")
    #is_prescreened = True

    is_n50_ok = (
        n50 == asm_len  # single-sequence
        or n50 > 20000
        or ((gcf or gx_kingdom == "prst") and n50 > 10000)
        or blast_div == "viruses"
    )

    status = (
        "not-in-scope"         if not ftp_path
        else "ok-acceptlist"   if is_accept
        else "denylist"        if is_reject
        else "not-in-scope"    if weight < 2 or gx_div == "unkn:unknown"
        else "too-large"       if asm_len > 4.3e9
        else "low-n50"         if not is_n50_ok
        else "not-prescreened" if not is_prescreened
        else "ok"
    )

    num_bioprojects = tax2num_bioprojects.get(tax_id, 0)

    return Assembly(
        tax_id, asm_len, n50, species, orgname, gca, gcf, blast_div, gx_div,
        asm_name, asm_id, rs_category, excl_from_refseq, weight, num_bioprojects, sp_tax_id, tax_id, status, ftp_path,
    )
    # fmt: on


#############################################################################
def load_assemblies() -> List[Assembly]:
    assemblies = []
    fetch_assemblies_esummaries(out_filename="assemblies.esummary.xtract.tsv")
    with open("assemblies.esummary.xtract.tsv", "rt", encoding="utf8") as f:
        for line in f:
            line = line.rstrip("\n")

            if not line or line[0] == "#":
                continue
            assemblies.append(line2asm(line))

    assert len(assemblies) > 110000
    return assemblies


#############################################################################
# GP-33261
def filter_close_tax_neighbors(assemblies: List[Assembly]):  # -> sequence-of-Assembly

    # We have some very well represented blast-divs, so for those
    # we want to prune more aggressively using lower thresholds.
    # fmt: off
    tax_proximity_thresholds = {
        "sharks and rays"   : 0.25,
        "bony fishes"       : 0.25,  # this selects 25 fishes
        "eudicots"          : 0.35,
        "rodents"           : 0.38,  # GP-33335
        "crustaceans"       : 0.40,
        "birds"             : 0.45,  # 0.4 selects 6 birds; 0.45 selects 25
        "ascomycete fungi"  : 0.40,
        "placentals"        : 0.40,
        "bats"              : 0.45,
        "turtles"           : 0.45,
        "primates"          : 0.55,  # GP-33335
    }
    # fmt: on
    proximity_thr = 0.55  # default
    assert all((div in blast_div2gx for div in tax_proximity_thresholds))

    neighbors = defaultdict(list)  # tax-id -> list of neighbor taxa.
    with gzip.open(args.tax_distances, "rt", encoding="ascii") as f:
        for line in f:
            if not line or line[0] == "#":
                continue

            cols = line.rstrip().split("\t")
            t1 = TaxId(int(cols[0]))
            t2 = TaxId(int(cols[1]))

            assert cols[2] in ("p")

            thr1 = tax_proximity_thresholds.get(str(tax2blast_div.get(t1)), proximity_thr)
            thr2 = tax_proximity_thresholds.get(str(tax2blast_div.get(t2)), proximity_thr)
            if t1 != t2 and float(cols[3]) > min(thr1, thr2):
                neighbors[t1].append(t2)
                neighbors[t2].append(t1)

    #########################################################################

    # taxa out of scope based on proxmity: tax-id -> already-selected close neighbor tax-id
    seen_taxa: Dict[TaxId, TaxId] = {}

    taxa_counts = Counter((a.tax_id for a in assemblies))
    assemblies.sort(reverse=True, key=lambda a: (
        a.weight,               # sort by our definition of get_weight()
        a.num_bioprojects,      # then by popularity based on number of bioprojects
        taxa_counts[a.tax_id],  # then by popularity based on count of assemblies for the taxon
        a.n50,                  # then by higher N50
        a.gca))                 # then by newer GCA-acc

    for a in assemblies:
        seen_tax_id = seen_taxa.get(a.tax_id) or seen_taxa.get(a.sp_tax_id)

        if not seen_tax_id and a.status == "not-prescreened":
            # Could use it (not yet seen), but needs not prescreened.
            # Change status to needs-prescreen.
            a = a._replace(status="needs-prescreen")

        elif seen_tax_id and a.status == "ok":
            # Already have this or close tax-id - change status to seen.
            a = a._replace(status="seen", nbr_tax_id=seen_tax_id)

        elif a.status.startswith("ok"):
            # Not-yet-seen, or on accept-list, or reference: update seen.
            seen_taxa[a.tax_id] = a.tax_id
            seen_taxa[a.sp_tax_id] = a.sp_tax_id
            for nbr_tax_id in neighbors[a.tax_id]:
                seen_taxa[nbr_tax_id] = a.tax_id

        yield a


#############################################################################
tot_len = 0.0
printed_header = False
counts_by_div: Dict[GxDiv, int] = defaultdict(int)

with open("assemblies.tsv", "wt", encoding="utf8") as f:
    for a in filter_close_tax_neighbors(load_assemblies()):
        if not printed_header:
            print("#", end="", file=f)
            print(*a._fields, sep="\t", file=f)
            printed_header = True

        if a.status.startswith("ok"):
            tot_len += a.len
            counts_by_div[a.gx_div] += 1

        print(*a, sep="\t", file=f)

eprint(f"Generated {args.out_dir}/assemblies.tsv")

#############################################################################

# fmt: off
# Verify expected assembly counts for top-20 divs
expected_asm_counts = {
    "virs:prokaryotic viruses"  : 14498,
    "virs:eukaryotic viruses"   : 7867,
    "prok:high GC Gram+"        : 3033,
    "prok:firmicutes"           : 2762,
    "prok:g-proteobacteria"     : 2544,
    "prok:a-proteobacteria"     : 2235,
    "prok:CFB group bacteria"   : 1798,
    "prok:bacteria"             : 1278,
    "fung:ascomycetes"          : 1224,
    "prok:b-proteobacteria"     : 937,
    "fung:basidiomycetes"       : 470,
    "arch:archaea"              : 453,
    "fung:budding yeasts"       : 432,
    "prok:d-proteobacteria"     : 262,
    "prok:proteobacteria"       : 200,
    "anml:insects"              : 193,
    "plnt:plants"               : 108,
    "prst:algae"                : 83,
    "anml:nematodes"            : 80,
}
# fmt: on

for gx_div, expected_count in expected_asm_counts.items():
    count = counts_by_div[gx_div]
    if count < expected_count * 0.9:
        eprint(
            f"Warning: count of {gx_div} assemblies decreased by more than 10%.",
            f"Expected: {expected_count}; Actual: {count}",
        )

eprint(
    "\nSummary of top-50 divs by total sequence: ", "(columns: blast-div, num-taxa, num-rs-asms, total-len, median-len)"
)

# $18 ~/ok/ - take selected.
# -g8 - group by blast-div, count tax-id, countunique GCA, sum len, median len
# sort by total length in blast div, descending
run(
    """ cat assemblies.tsv |
         awk -v FS='\\t' '$18 ~/ok/' |
         datamash -s -g8 count 1 countunique 7 sum 2 median 2 |
         sort -t$'\\t' -k4,4nr |
         head -n50 |
         column -t -s$'\\t'
    """
)

eprint("\nTotal selected length:", int(tot_len * 100 / 1e9) / 100, "Gbp")

# Sanity checks.
assert 680e9 < tot_len < 750e9
