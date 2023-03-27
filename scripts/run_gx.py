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

import sys
assert sys.version_info.major >= 3 and sys.version_info.minor >= 8, f"Python version: {sys.version_info}. Require python 3.8 or newer."

import argparse
import os
import signal
import re
import subprocess
import shlex
import shutil
import time
import json
import glob
import urllib.request

from typing import Tuple

from collections import defaultdict
from collections import namedtuple
from pathlib import Path

# Execute gx-pipeline:
# python3 run_gx.py
#   --tax-id=195
#   --fasta=/am/ftp-genomes/all/GCA/001/490/155/GCA_001490155.1_CCN120/GCA_001490155.1_CCN120_genomic.fna.gz
#   --div="e-proteobacteria"
#   --split-fasta=T
#   --out-dir=.

horizontal_line = "\n--------------------------------------------------------------------\n"


# ---------------------------------------------------------------------------
def eprint(*args):
    print(*args, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
g_verbose = False
g_use_color = bool(os.getenv("FCS_USE_COLOR"))


# Will print exception-text underlined, and in red if using colors
def excepthook(etype, value, tb):
    msg = "\033[4m" + ("\033[91m" if g_use_color else "") + str(value) + "\033[0m" + "\n"
    sys.__excepthook__(etype, etype(msg), tb)

sys.excepthook = excepthook


def debug_print(*args):
    if not g_verbose:
        pass
    elif g_use_color:
        eprint("\033[2m####### ", *args, "\033[0m") # dim mode
    else:
        eprint("#######", *args)


def with_this_py(args):
    # args[0] may be a bazel executable (without .py), or a python script (with .py)
    # If the latter, execute it with same python3 binary as this script, 
    # in case the user invoked it with a specific python path, overriding shebang.
    if os.path.exists(args[0] + ".py"):
        args[0] += ".py"

    if args[0].endswith(".py"):
        args = [sys.executable] + args

    return args


# ---------------------------------------------------------------------------
class ProcessPipeline:
    """
    Emulate shell pipelines and process-substitutions in pure python and anonymous pipes.

    Each ProcessPipeline represents a chain of `cmd | cmd2 | ... | cmdN [ >out_file ]`
    Such chains are composable to form arbitrary dataflows.

    Consider the following example Bash pipeline:

        diff <(zcat -f file1.txt | grep -v foo | sort)
             <(zcat -f file2.txt | sort)    |
        tee >(gzip -c >file1_file2.diff.gz) |
        wc > wc.txt

    It would be implemented as follows:

    with (
        ProcessPipeline() as p1,  # arg-1 of diff
        ProcessPipeline() as p2,  # arg-2 of diff
        ProcessPipeline() as p3,  # arg of tee
        ProcessPipeline() as p0,  # diff | tee | wc > wc.txt
    ):
        p1.add("zcat -f file1.txt")
        p1.add("grep -v foo", ok_returncodes=(0,1)) # grep returns 1 if no matches
        p1.add("sort")

        p2.add("zcat -f file2.txt")
        p2.add("sort")

        p3.add("gzip -c", stdin=subprocess.PIPE, out_filename="file1_file2.diff.gz")

        # diff returns 0 if no diffs, 1 if there are diffs, >1 if error.
        p0.add(f"diff /dev/fd/{p1.stdout_fd} /dev/fd/{p2.stdout_fd}", ok_returncodes=(0,1))
        p0.add(f"tee /dev/fd/{p3.stdin_fd}")
        p0.add("wc", out_filename="wc.txt")
    """

    # p              : process' handle returned by subprocess.Popen.
    # cmd            : command used to start the process.
    # ok_returncodes : tuple of allowed ret-code values.
    # fds            : tuple of fds to close on exit.
    Process = namedtuple("Process", "p cmd ok_returncodes fds")  # pylint: disable=C0103

    def __init__(self):
        self.processes = []
        self.out_file = None

    @property
    def stdin_fd(self):
        return self.processes[0].p.stdin.fileno()

    @property
    def stdout_fd(self):
        return self.processes[-1].p.stdout.fileno()

    # -----------------------------------------------------------------------
    # Add a pipeline stage. Connect last stage's stdout to this stage's stdin.
    def add(self, cmd, ok_returncodes=(0,), out_filename=None, stdin=None) -> None:

        debug_print("Starting process", cmd)

        assert not (
            self.processes and stdin
        ), f"stdin can only be specified for the first stage. {self.processes}, {stdin}"

        if out_filename:
            assert not self.out_file, "Already have terminal stage with output-file."
            self.out_file = open(out_filename, "wb")  # pylint: disable=R1732

        fds = tuple(int(fd) for fd in re.findall(r"\W/dev/fd/(\d+)", str(cmd)))

        p = subprocess.Popen(  # pylint: disable=R1732
            shlex.split(cmd) if isinstance(cmd, str) else cmd,
            stdout=self.out_file if self.out_file else subprocess.PIPE,
            stdin=(stdin if stdin else self.processes[-1].p.stdout if self.processes else None),
            stderr=sys.stderr,
            pass_fds=fds,
        )

        self.processes.append(self.Process(p, str(cmd), ok_returncodes, fds))

        if len(self.processes) > 1:
            # For SIGPIPE backpropagation.
            # See https://docs.python.org/3/library/subprocess.html#replacing-shell-pipeline
            self.processes[-2].p.stdout.close()

        return self

    # -----------------------------------------------------------------------
    @staticmethod
    def cleanup(processes) -> int:  # Return number of processes exited with error
        num_errors = 0

        for p in processes:
            debug_print("Cleaning up process", p.cmd)
            p.p.terminate()

            # closing fd signals to the process on the other end
            # of anonymous pipe that the pipe is closed, unblocking it.
            for fd in p.fds:
                os.close(fd)

            if p.p.returncode is not None and p.p.returncode not in p.ok_returncodes:
                eprint(f"Error: Process failed with retcode {p.p.returncode}: {p.cmd})")
                num_errors += 1

        return num_errors

    # -----------------------------------------------------------------------
    # NB: with processes connected by named pipes there's a following issue to
    # keep in mind: if we're waiting on a process (or a pipeline) that's connected
    # to another process via a named-pipe, and the process on the other end of
    # the pipe errors-out, our wait blocks indefinitely. So the wait-loop must
    # monitor ALL processes involved, not just this particular pipeline.
    #
    # With anonymous pipes this is not an issue becaus the OS closes the pipe
    # when the other process errors-out, so we can wait only on pipe's processes.
    def wait(self) -> None:
        while self.processes:
            time.sleep(0.01)
            cleanup_targets = [p for p in self.processes if p.p.poll() is not None]
            num_errors = ProcessPipeline.cleanup(cleanup_targets)
            self.processes = [p for p in self.processes if p not in cleanup_targets]
            assert num_errors == 0, "Had errors."

    # -----------------------------------------------------------------------
    def __enter__(self):
        return self

    # -----------------------------------------------------------------------
    def __exit__(self, exc_type, exc_value, traceback):
        try:
            if not traceback:
                self.wait()
        finally:
            ProcessPipeline.cleanup(self.processes)
            self.processes = []
            if self.out_file:
                self.out_file.close()
                self.out_file = None


# ---------------------------------------------------------------------------
def get_json(url):
    with urllib.request.urlopen(url, timeout=10) as r:
        assert r.status == 200, f"Error: {url} - status {r.status}."
        return json.loads(r.read().decode(r.headers.get_content_charset()))


# ---------------------------------------------------------------------------
eutils_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get_assembly_esummary(gc_acc):
    j = get_json(f"{eutils_url}/esearch.fcgi?db=assembly&format=json&term={gc_acc}%5BAssembly%20Accession%5D")
    # NB: may be more than one - GP-33342
    assert len(j["esearchresult"]["idlist"]) >= 1, f"Did not find assembly {gc_acc}"
    asm_uid = int(j["esearchresult"]["idlist"][0])

    j = get_json(f"{eutils_url}/esummary.fcgi?db=assembly&format=json&id={asm_uid}")
    return j["result"][str(asm_uid)]


# ---------------------------------------------------------------------------
def get_tax_id_from_species(species):
    if not species:
        return None

    sp = urllib.parse.quote(species)
    j = get_json(f"{eutils_url}/esearch.fcgi?db=taxonomy&format=json&term={sp}%5Bporgn%5D")
    assert len(j["esearchresult"]["idlist"]) == 1, f"Did find unique tax-id for '{species}'"
    return int(j["esearchresult"]["idlist"][0])


# ---------------------------------------------------------------------------
def redirect_stdout_stderr_to_file(args) -> None:
    pid = os.getpid()
    filename = f"{args.out_dir}/tmp_{pid}.summary.txt"
    sys.stdout = open(filename, 'w')
    sys.stderr = sys.stdout
    return filename


# ---------------------------------------------------------------------------
def guess_assembly_path(gc_acc: str, gc_genomes_root_dir: str):
    ps = [p for p in Path("/").glob(gc_genomes_root_dir.strip("/") + f"*/all/{gc_acc[0:3]}/{gc_acc[4:7]}/{gc_acc[7:10]}/{gc_acc[10:13]}/{gc_acc}_*/{gc_acc}_*_genomic.fna.gz") if "_cds_from_" not in str(p) and "_rna_from_" not in str(p)]
    return str(ps[0]) if len(ps) == 1 else None


# ---------------------------------------------------------------------------
# Get (tax_id, detailed-basename, ftp-fasta-path) from Eutils Assembly summary.
# e.g. (
#   8407,
#   "GCF_905171775.1__8407__Rana_temporaria__common_frog",
#   "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/330/725/GCF_003330725.1_ASM333072v1/GCF_003330725.1_ASM333072v1_genomic.fna.gz" pylint: disable=C0301
# )
def get_taxid__basename__ftp_path(gc_acc: str, gc_genomes_root_dir: str) -> tuple:
    asm_info = get_assembly_esummary(gc_acc)
    tax_id = int(asm_info["taxid"])

    org = asm_info["organism"]  # ......... e.g. Acidisarcina polymorpha (bacteria)
    org = re.sub(r"[^\w.-]", "_", org)  # . e.g. Acidisarcina_polymorpha__bacteria_
    org = re.sub(r"_+$", "", org)  # ...... e.g. Acidisarcina_polymorpha__bacteria

    # NB: Do not use asm_info["ftppath_assembly_rpt"] - see GP-33326
    ftp_path = asm_info["ftppath_genbank" if gc_acc.startswith("GCA_") else "ftppath_refseq"]
    assert ftp_path, "No FTP path for this assembly."
    assert not ftp_path.endswith("/")
    ftp_path += "/" + ftp_path.split("/")[-1] + "_genomic.fna.gz"

    # GP-32574
    prefix = "ftp://ftp.ncbi.nlm.nih.gov/genomes"
    assert ftp_path.startswith(prefix), (ftp_path, prefix)

    if gc_genomes_root_dir:
        # Convert ftp-path to filesystem-path
        assert Path(gc_genomes_root_dir).is_dir(), f"Can't access --gc-genomes-root-dir={gc_genomes_root_dir}"
        paths = [gc_genomes_root_dir.rstrip("/") + sfx + ftp_path[len(prefix) :] for sfx in ["", "2", "3", "4", "5"]]
        ftp_path = next((p for p in paths if Path(p).is_file()), None)
        if not ftp_path:
            # Don't assert here; just warn - will fall-back on guess_assembly_path in the caller code.
            eprint(f"Can't resolve fasta path from FTP path {paths}")

    return (tax_id, f"{gc_acc}__{tax_id}__{org}", ftp_path)


# ---------------------------------------------------------------------------
# Return first line from file matching pattern.
#
# NB: The pythonic way is about 10x slower in comparison -
# ~1.7 seconds instead of ~0.25 with gzip or ~0.12 with minigzip, which adds up for small assemblies.
# with gzip.open(args.gx_db.replace(".gxi", ".blast_div.tsv.gz"), "rt", encoding="ascii") as f:
#     ret = next((line.rstrip().split("\t")[1] for line in f if line.startswith(str(args.tax_id) + "\t")), None)
# """
def grep(file_path: str, pattern: str):
    if not Path(file_path).is_file():
        return None

    with ProcessPipeline() as p:
        p.add(["zcat", "-f", file_path])
        p.add(["grep", "-E", pattern], ok_returncodes=(0, 1))
        return p.processes[-1].p.communicate()[0].decode("ascii")


# ---------------------------------------------------------------------------
# Look-up blast-div for tax-id from *blast_di.tsv.gz, with eutils Taxonomy as fallback.
def get_blast_div(args) -> str:
    assert args.tax_id
    assert args.gx_db

    line = grep(args.gx_db.replace(".gxi", ".blast_div.tsv.gz"), f"^{args.tax_id}\t")
    div = line.split("\t")[1].rstrip() if line else None

    # Fallback on eutils:
    def get_div_from_taxonomy():
        debug_print("Getting div from Taxonomy.")
        j = get_json(f"{eutils_url}/esummary.fcgi?db=taxonomy&format=json&id={args.tax_id}")
        try:
            return j["result"][str(args.tax_id)]["division"]
        except:
            eprint(f"\nError: Could not look-up tax-id {args.tax_id} from NCBI Taxonomy.\n\n")
            raise

    return div or get_div_from_taxonomy()


# -----------------------------------------------------------------------
def check_preconditions(args) -> None:

    if not Path(f"{args.gx_db}").is_file():
        eprint(f"GX-database {args.gx_db} is not mounted on this host ({os.uname().nodename}).\n")
        instructions_md = f"{args.bin_dir}/gxdb_cp_instructions_ncbi.md"
        if args.gx_db.startswith("/dev/shm/gxdb") and Path(instructions_md).is_file():
            subprocess.run(["cat", instructions_md], check=True)
            sys.exit(2)

    for file in ["gx", "classify_taxonomy", "action_report", "blast_names_mapping.tsv"]:
        ok = any((Path(f"{args.bin_dir}/{file}{sfx}").is_file() for sfx in ("", ".py")))
        assert ok, f"File {args.bin_dir}/{file} is not accessible. Might need to specify --bin-dir explicitly."


# ---------------------------------------------------------------------------
# Set fasta, tax_id, out_basename, div args, if not specified.
# Also set out_taxonomy_rpt and gzip_c args.
# Verify that files are accessible.
def fill_missing_args(args) -> None:

    # minigzip is a faster (not quite drop-in) replacement for gzip - use it if available.
    # NB: minigzip -dc /path/to/file does not work - must cat | minigzip ... instead
    args.gzip_c = "minigzip -c" if shutil.which("minigzip") else "gzip -c"

    # GP-33682
    if not args.fasta and args.gc_genomes_root_dir and args.gc_acc:
        args.fasta = guess_assembly_path(args.gc_acc, args.gc_genomes_root_dir)

    # -----------------------------------------------------------------------
    # Fill-in missing args computable from gc-acc.
    (gc_tax_id, gc_basename, gc_ftp_fasta_path) = (
        get_taxid__basename__ftp_path(args.gc_acc, args.gc_genomes_root_dir)
        if args.gc_acc and not (args.out_basename and args.tax_id and args.fasta)
        else (None, None, None)
    )

    species_tax_id = get_tax_id_from_species(args.species)

    args.tax_id = args.tax_id or gc_tax_id or species_tax_id
    assert args.tax_id, "--species or --tax-id or --gc-acc must be provided."

    # if derived tax-id by multiple means, verify that values are consistent
    assert not gc_tax_id      or gc_tax_id      == args.tax_id, (gc_tax_id, args.tax_id)
    assert not species_tax_id or species_tax_id == args.tax_id, (species_tax_id, args.tax_id)

    args.fasta = args.fasta or gc_ftp_fasta_path
    assert args.fasta and (
        args.fasta.startswith("ftp://") or Path(args.fasta).exists()  # not is_file() - can be /dev/stdin
    ), "--fasta or --gc-genomes-root-dir and --gc-acc must be provided."

    args.out_basename = args.out_basename or gc_basename

    eprint(horizontal_line)
    eprint("tax-id    :", args.tax_id)
    eprint("fasta     :", args.fasta)
    if Path(args.fasta).is_file() and not args.fasta.endswith(".mft"):
        eprint("size      :", int(os.path.getsize(args.fasta) * 100 / (1024 * 1024)) / 100, "MiB")
    eprint("split-fa  :", args.split_fasta)

    # -----------------------------------------------------------------------
    args.div = args.div or get_blast_div(args)
    assert args.div, "Could not look-up blast-div for the tax-id. Must be specified wih --div"

    # -----------------------------------------------------------------------
    # Convert blast-div to gx-div (unless already gx-div)
    if len(args.div) < 5 or args.div[4] != ":":
        eprint("BLAST-div :", args.div)
        gx_div = None
        with open(f"{args.bin_dir}/blast_names_mapping.tsv", "rt", encoding="ascii") as f:
            gx_div = next((line.rstrip().split("\t")[1] for line in f if line.startswith(args.div + "\t")), None)

        if not gx_div:
            eprint(f">>>>>> Warning: Taxonomy div '{args.div}' is not known to GX. using 'unknown'.")
            eprint(">>>>>> NB: GX-div can be specified directly, e.g. --div=anml:molluscs")
            gx_div = "unkn:unknown"

        args.div = gx_div

    assert args.div[4] == ":", "Expected gx-div with kingdom-prefix."
    eprint("gx-div    :", args.div)

    # -----------------------------------------------------------------------
    # GP-33305
    # If allow_same_species is unspecified, set to True if gc-acc is unspecified
    # or the assembly-accesison-no-version is not in the index scope, else False.
    if args.allow_same_species is None:
        args.allow_same_species = not (args.gc_acc and grep(
            args.gx_db.replace(".gxi", ".assemblies.tsv"),
            "\t" + args.gc_acc.split(".")[0].replace("GCA_", "GC[AF]_").replace("GCF_", "GC[AF]_")
        ))
    eprint("w/same-tax:", args.allow_same_species)

    # -----------------------------------------------------------------------
    # Construct out_basename, if not provided, and prefix out_dir.
    assert args.fasta != "/dev/stdin" or args.out_basename, "--out-basename must be specified."
    args.out_basename = args.out_basename or re.sub(r"\.\w+$", f".{args.tax_id}", Path(args.fasta).name)
    args.out_basename = f"{args.out_dir}/{args.out_basename}"

    args.out_taxonomy_rpt = f"{args.out_basename}.taxonomy.rpt"
    eprint("bin-dir   :", args.bin_dir)
    eprint("gx-db     :", args.gx_db)

    gx_help = subprocess.run([f"{args.bin_dir}/gx", "--help"], check=True, capture_output=True, encoding="ascii")
    gx_version = [line[6:] for line in gx_help.stdout.split("\n") if line.startswith("build:")][0].strip()
    eprint("gx-ver    :", gx_version)

    eprint("output    :", args.out_taxonomy_rpt)

    eprint(horizontal_line)

    debug_print("args:", args, "\n")


# ---------------------------------------------------------------------------
# Execute the GX pipeline.
# Equivalent of:
#
# zcat -f $fasta |
#    /path/to/gx split-fasta |
#    pv -Wbrat --interval 0.5 |
#    /path/to/gx align --db=/path/to/gxdb/all.gxi --repeats-basis-fa=<(zcat -f $fasta) |
#    tee >(gzip -c > $out_basename.hits.tsv.gz) |
#    awk -v FS='\t' -v tax_id=$tax_id '($3 != tax_id)' |
#    /path/to/gx taxify --db=/path/to/gxdb/all.gxi -o $out_basename.taxonomy.rpt.tmp
def run_gx_pipeline(args) -> None:
    def add_zcat_fasta(p):
        p.add(["gzip", "-cdf", args.fasta])
        if args.fasta.endswith(".mft"):
            p.add(["grep", "-Ev", "^(#|$)"])
            p.add(["xargs", "-n1", "gzip", "-cdf"])

    def run(p_zcat_fasta, p_save_hits, p_main):
        Path(args.out_dir).mkdir(parents=True, exist_ok=True)

        # Download fasta, as appropraite
        if args.fasta.startswith("ftp://"):
            local_fasta = args.out_dir + "/" + args.fasta.split("/")[-1]
            if not Path(local_fasta).exists():
                eprint(f"Downloading fasta to {local_fasta}...")
                urllib.request.urlretrieve(args.fasta, local_fasta + ".download")
                os.rename(local_fasta + ".download", local_fasta)
            args.fasta = local_fasta

        gx_bin = f"{args.bin_dir}/gx"

        add_zcat_fasta(p_main)

        if args.split_fasta:
            p_main.add([gx_bin, "split-fasta"])

        if shutil.which("pv"):
            # For .mft files size is unknown, so will use --size=0 for simplified bar and no ETA
            factor = (3.2 if args.fasta.endswith(".gz") else 0 if args.fasta.endswith(".mft") else 1)
            est_bytes = int(os.path.getsize(args.fasta) * factor)
            p_main.add(["pv", "-Wbratpe", "--interval=0.5", f"--size={est_bytes}"])

        # For large genomes (not prok/virus/synthetic) will mask transposon repeats (unless explicitly specified).
        with_repeats = (
            args.mask_transposons if args.mask_transposons is not None
            else not any((args.div.startswith(s) for s in ["synt", "prok", "arch", "virs", "unkn"]))  # NB: unkn is for metagenomes
        )

        assert not (with_repeats and args.fasta == "/dev/stdin"), "Need multiple passes over the input, so --fasta must be a file"

        if with_repeats:
            add_zcat_fasta(p_zcat_fasta)

        # NB: time and /usr/bin/time are different things.
        p_main.add(
            []
            + ([shutil.which("time"), "-v"] if args.debug and Path(shutil.which("time")).is_file() else [])
            + (["nice", "-n19"] if shutil.which("nice") else [])
            + [gx_bin, "align", f"--gx-db={args.gx_db}"]
            + ([f"--repeats-basis-fa=/dev/fd/{p_zcat_fasta.stdout_fd}"] if with_repeats else [])
        )

        # As-if tee >(gzip -c >{args.out_basename}.hits.tsv.gz)
        if args.save_hits:
            p_save_hits.add(args.gzip_c, stdin=subprocess.PIPE, out_filename=f"{args.out_basename}.hits.tsv.gz")
            p_main.add(["tee", f"/dev/fd/{p_save_hits.stdin_fd}"])

        if not args.allow_same_species:
            p_main.add(["awk", "-v", "FS=\t", f"(NR>2 && NF != 8){{ exit 111; }} ($3 != {args.tax_id})"])

        p_main.add([
            gx_bin,
            "taxify", 
            f"--gx-db={args.gx_db}", 
            f"--output={args.out_taxonomy_rpt}.tmp",
            f"--asserted-div={args.div}",
        ] + [f"--db-exclude-locs={path}" for path in (os.getenv("GX_EXCLUDE_LOCS", f"{args.bin_dir}/db_exclude.locs.tsv"),) if os.path.exists(path)])  # GP-34552

    # NB: nested-with instead of multi-with to support python3.8
    with ProcessPipeline() as p_zcat_fasta:
        with ProcessPipeline() as p_save_hits:
            with ProcessPipeline() as p_main:
                run(p_zcat_fasta, p_save_hits, p_main)


# ---------------------------------------------------------------------------
# Run classify_taxonomy and action_report
def run_classify_taxonomy_and_action_report(args) -> None:
    def run(cmd, out_filename):
        with open(out_filename, "wb") as out_file:
            subprocess.run(cmd, stdout=out_file, check=True, stderr=sys.stderr)

    run(
        with_this_py([f"{args.bin_dir}/classify_taxonomy", f"--in={args.out_taxonomy_rpt}.tmp"]),
        args.out_taxonomy_rpt,
    )
    os.remove(f"{args.out_taxonomy_rpt}.tmp")

    if args.action_report:
        run(
            with_this_py([f"{args.bin_dir}/action_report", f"--in={args.out_taxonomy_rpt}"]),
            f"{args.out_basename}.fcs_gx_report.txt",
        )


# ---------------------------------------------------------------------------
def print_summary(
    # fmt: off
    args,
    file     : str,  # (.taxonomy.rpt or .fcs_gx_report.txt)
    title    : str,  # "contamination" or "action"
    grp_col  : int,  # column on which to group on
    coords   : Tuple[int, ...],  # (len) col, or (start,stop) cols
    predicate_fn,
    # fmt: on
) -> None:

    header = f"{file} {title} summary:"
    eprint(f"{horizontal_line}\n{header}\n" + (len(header) * "-"))

    lens = defaultdict(int)
    counts = defaultdict(int)

    with open(f"{args.out_basename}.{file}", "r", encoding="ascii") as f:
        for line in f:
            if not line or line[0] == "#":
                continue

            cols = line.rstrip().split("\t")
            assert len(cols) in (8, 34)  # fcs_genome.txt or taxonomy.rpt

            if not predicate_fn(cols):
                continue

            counts[cols[grp_col]] += 1
            lens[cols[grp_col]] += (
                int(cols[coords[0]]) if len(coords) == 1 else (int(cols[coords[1]]) - int(cols[coords[0]]) + 1)
            )

    arr = [(div, counts[div], len) for (div, len) in lens.items()]
    arr.sort(key=lambda x: x[2], reverse=True)
    arr.insert(0, ("-----", "-----", "----------"))
    arr.insert(0, ("TOTAL", sum(counts.values()), sum(lens.values())))
    arr.insert(0, ("", "-----", "----------"))
    arr.insert(0, ("", "seqs", "bases"))
    for row in arr:
        eprint("{: <30} {: >5} {: >10}".format(*row))  # pylint: disable=C0209


# ---------------------------------------------------------------------------
def parse_args():

    # --gc-acc will only be available if is_NCBI
    is_ncbi = bool(os.getenv("NCBI", ""))

    parser = argparse.ArgumentParser(
        description="Execute GX-pipeline (gx-split-fasta|gx-align|gx-taxify), classify-taxonomy, and genome-report steps.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # to show defaults
        prog="screen genome"  # GP-34867 (run_gx will be aliased as "screen genome" in the fcs.py runner script)
    )

    parser._action_groups.pop()  # pylint: disable=W0212
    req = parser.add_argument_group("required arguments (if --gc-acc not specified).")
    opt = parser.add_argument_group("optional inputs")
    out = parser.add_argument_group("outputs")
    other = parser.add_argument_group("other")

    def str2bool(x):
        ts = ("yes", "true", "t", "y", "1")
        fs = ("no", "false", "f", "n", "0")
        x = x.lower()

        if x in ("none", ""):
            return None

        assert x in ts or x in fs, f"{x}: Boolean value expected."
        return x in ts

    # -----------------------------------------------------------------------
    # Required args.

    req.add_argument(
        "--fasta",
        metavar="FASTA_FILE",
        required=(not is_ncbi),
        help="Input FASTA file; may be gzipped; may be a manifest-file (with .mft extension), or ftp-path.",
    )

    req.add_argument(
        "--tax-id",
        type=int,
        required=(not is_ncbi),
        help="NCBI Taxonomy identifier.",
    )

    req.add_argument(
        "--species",
        type=str,
        required=False,
        help="Alternatively or in addition to --tax-id, can specify species binomial name.",
    )

    # -----------------------------------------------------------------------
    # Optional args.

    if is_ncbi:
        opt.add_argument(
            "--gc-acc",
            help="NCBI Assembly accession.",
        )

        opt.add_argument(
            "--gc-genomes-root-dir",
            default=os.getenv("GC_GENOMES_ROOT_DIR"),
            help="GenColl assemblies FTP root dir.",
        )

    opt.add_argument(
        "--split-fasta",
        type=str2bool,
        metavar="BOOL",
        default=True,
        help="Split fasta sequences on N-runs.",
    )

    opt.add_argument(
        "--div",  # Can also be gx-div, e.g. anml:primates
        metavar="DIV",
        help="BLAST-div of the tax-id, from 'NCBI BLAST name' on taxon Info page.",
    )

    opt.add_argument(
        "--gx-db",
        default=os.getenv("GX_DB_DIR", "/dev/shm/gxdb/"),
        help="Path to gx-db, including the filename-prefix.",
    )

    opt.add_argument(
        "--mask-transposons",
        type=str2bool,
        metavar="BOOL",
        help="Whether to mask transposons in the input. If not specified, will mask for euks only.",
    )

    def get_default_bin_dir():
        # NB, prefer $_ env-var because in bazel-binary sys.argv[0] is
        # is a path into /tmp/Bazel.runfiles/.../run_gx, not the path of the
        # bazel-binary itself, whereas _ env-var is.
        # https://unix.stackexchange.com/questions/292991
        #
        # However, if the script is called as `python run_gx.py`, then $_
        # is the python interpreter's path, not that of the script.
        #
        # Additionally, the path could be a symlink, and the bin-dir
        # is that of the symlink's target, or that of the symlink.

        paths = [os.getenv("_"), sys.argv[0]]
        paths += [os.readlink(p) for p in paths if p and Path(p).is_symlink()]
        dirs = [os.path.dirname(p) for p in paths if p]  # NB: p can be None
        return next((d for d in dirs if Path(d + "/gx").is_file()), ".")

    opt.add_argument(
        "--bin-dir",
        default=os.getenv("GX_BIN_DIR", default=get_default_bin_dir()),
        help="Path to executables.",
    )

    opt.add_argument(
        "--allow-same-species",
        type=str2bool,
        metavar="BOOL",
        default=None,
        help="Whether to use same-species hits as evidence.",
    )
    
    # -----------------------------------------------------------------------
    # Outputs.

    out.add_argument(
        "--out-basename",
        metavar="NAME",
        help="Output files' basename. Default is {fasta-basename}.{tax-id}.",
    )

    out.add_argument(
        "--out-dir",
        default=".",
        help="Output directory.",
    )

    out.add_argument(
        "--action-report",
        type=str2bool,
        metavar="BOOL",
        default=True,
        help="Generate *.fcs_action_report.txt output.",
    )

    out.add_argument(
        "--save-hits",
        type=str2bool,
        metavar="BOOL",
        default=False,
        help="Save intermediate alignments into {out-basename}.hits.tsv.gz.",
    )

    out.add_argument(
        "--generate-logfile",
        type=str2bool,
        metavar="BOOL",
        default=False,
        help="Redirect stdout and stderr to {out-basename}.summary.txt.",
    )

    # -----------------------------------------------------------------------

    other.add_argument(
        "--debug",
        action="store_true",
        help="Enable verbose mode.",
    )

    # print help and error-out if no args are supplied
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)

    args = parser.parse_args()

    # These options were added conditionally,
    # so set them to None so we don't have to check hasattr later.
    if not is_ncbi:
        args.gc_acc = None
        args.gc_genomes_root_dir = None

    elif not args.gc_genomes_root_dir and args.fasta and args.fasta.startswith("ftp://"):
        eprint(
            "\n-----------------------------------------",
            "\n NOTE: to avoid downloading fasta from FTP,",
            "\n please specify --gc-genomes-root-dir=...",
            "\n or set GC_GENOMES_ROOT_DIR env-var.",
            "\n-----------------------------------------",
            "\n")

    return args


# -----------------------------------------------------------------------
def main() -> None:

    # To propagate signals to child processes:
    # https://stackoverflow.com/questions/67823770/how-to-propagate-sigterm-to-children-created-via-subprocess
    def signal_handler(signo, _):
        time.sleep(2)  # give the child processes a moment to wrap-up
        raise SystemExit(f"Signal: {signo} ({signal.strsignal(signo)})")

    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGQUIT):
        signal.signal(sig, signal_handler)

    args = parse_args()

    global g_verbose
    g_verbose = args.debug

    if args.generate_logfile:
        log_file_name = redirect_stdout_stderr_to_file(args)

        # NB1: we are redirecting to a temp-file first because 
        # basename is computed in fill_missing_args - it is not available here yet.
        # Will rename the file in the end.
        # 
        # NB2: The correct way to redirect is described below:
        # https://eli.thegreenplace.net/2015/redirecting-all-kinds-of-stdout-in-python/
        # 
        # The simple way should be enough for our purposes.

    if not args.gx_db.endswith(".gxi"):
        # Can be specified uniquely in any of the following ways:
        # /dev/shm/gxdb/all.gxi /dev/shm/gxdb/all /dev/shm/gxdb/ /dev/shm/gxdb
        # NB: this must be done before check_preconditions call below
        paths = glob.glob(args.gx_db + ("/" if os.path.isdir(args.gx_db) else "") + "*.gxi")
        assert len(paths) == 1 or os.path.exists(args.gx_db), f"Cannot access {args.gx_db}"
        assert len(paths) == 1, f"Cannot resolve path to *.gxi file from {args.gx_db}: {paths}"
        args.gx_db = paths[0]

    check_preconditions(args)
    fill_missing_args(args)
    run_gx_pipeline(args)
    run_classify_taxonomy_and_action_report(args)

    if args.action_report:
        print_summary(args, "fcs_gx_report.txt", "contamination", 5, (1, 2), lambda cols: cols[6] != "REVIEW" and cols[6] != "REVIEW_RARE")
        print_summary(args, "fcs_gx_report.txt", "action", 4, (1, 2), lambda cols: True)
    else:
        print_summary(args, "taxonomy.rpt", "contamination", 32, (1,), lambda cols: "contaminant" in cols[31])
    
    eprint(horizontal_line)

    # now that we have the the basename available, we can rename the log file appropriately
    if args.generate_logfile:
        new_log_file_name =  f"{args.out_basename}.summary.txt"
        os.rename(log_file_name, new_log_file_name)
        

# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
