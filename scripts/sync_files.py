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
# pylint: disable=W0622,C0301,C0116,C0103,C0114,C0115,W0511,R0912,R0914,R0915,R0903,W0702,W0603

import sys
import os
import argparse
import urllib.request
import subprocess
import json
import shutil
import ssl
import atexit
import time
import xml.etree.ElementTree as ET
import signal

from pathlib import Path
from collections import namedtuple
from datetime import datetime
from datetime import timedelta
#from multiprocessing import Pool
import concurrent.futures
import hashlib


# ---------------------------------------------------------------------------
def parse_args():
    is_ncbi = bool(os.getenv("NCBI", ""))

    # Will expose make-manifest command only if running within NCBI. GP-34867
    make_manifest_help_str = "make-manifest: Create the manifest-file from directory contents." if is_ncbi else ""

    parser = argparse.ArgumentParser(
        description=f"""File transfer utility.

    get: Make the content of DIR to be consistent with the input manifest.
         (transfer the files that are missing or differ; skip the files that are up-to-date).
         Note: the output directory is temporarily renamed to DIR.in_progress until done.

    check: Show which files in the manifest are missing or differ.

    {make_manifest_help_str}""",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="db"  # GP-34867 (sync_files.py will be aliased as "db get" "db check" etc)
    )

    parser._action_groups.pop()  # pylint: disable=W0212
    req = parser.add_argument_group('required arguments')

    req.add_argument(
        "command",
        choices=(["make-manifest"] if is_ncbi else[]) + ["get", "check"]
    )

    req.add_argument(
        "--mft",
        required=True,
        help="URL or rclone-remote:path or local-path to JSON-manifest or metalink (targets are relative to manifest's directory).",
    )

    req.add_argument(
        "--dir",
        required=True,
        default="",
        help="Output directory (or input in make-manifest mode)",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
class Color:
    use_color = bool(os.getenv("FCS_USE_COLOR"))
    (RESET, BOLD, DIM, UNDERLN, RED) = (
        ('\033[0m', '\033[1m', '\033[2m', '\033[4m', '\033[91m') if use_color 
        else ('', '', '', '', '')
    )

# ---------------------------------------------------------------------------
def eprint(*args, **kwargs):
    color_ = kwargs.get("color", Color.DIM)
    kwargs.pop("color", None)
    print(color_, file=sys.stderr, end="")
    print(*args, Color.RESET, file=sys.stderr, flush=True, **kwargs)


# ---------------------------------------------------------------------------
# Represent number of bytes as human-readable string with IEC-prefix: KiB, MiB, etc.
def as_readable(num) -> str:
    s = "PTGMK "
    while len(s) > 1 and abs(num) > 1024:
        (num, s) = (num / 1024, s[:-1])
    return f"{num}B" if s[-1] == " " else f"{int(num * 100 + 0.5) / 100}{s[-1]}iB"

assert as_readable(-123456) == "-120.55KiB"


# ---------------------------------------------------------------------------
def from_cmd(*args) -> str:
    result = subprocess.run(args, check=False, capture_output=True, encoding="ascii")
    # check=False - will check returncode ourselves and report command's stderr if it failed.
    if result.returncode != 0:
        eprint(result.stderr, color=Color.RED)
        assert False, (f"Command failed with retcode {result.returncode}:", args)
    return result.stdout


# ---------------------------------------------------------------------------
def get_content_from_ordinary_url(url) -> str:
    # ordinary meaning not gs://... s3://... file-path, etc.

    # Do we need this cert-workaround business?
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    with urllib.request.urlopen(url, timeout=10, context=ctx) as r:
        # NB: status can be either unset or None if url is ftp.
        assert not hasattr(r, "status") or r.status is None or r.status == 200, f"Error: {url} - status {r.status}."
        return r.read().decode(r.headers.get_content_charset() or "ascii")


# ---------------------------------------------------------------------------
def get_content(url) -> str:
    url_str = str(url)
    try :
        return (
                 from_cmd("gsutil", "cat", url_str)     if url_str.startswith("gs://")
            else get_content_from_ordinary_url(url_str) if "://" in url_str
            else from_cmd("rclone", "cat", url_str)     if ":/" in url_str
            else from_cmd("cat", url_str)
        ).strip()

    except:
        eprint(f"Cannot open {url}", color=Color.RED)
        raise


# ---------------------------------------------------------------------------
def remove_if_exists(path):
    if path and os.path.exists(path):
        eprint(f"\nRemoving {path}.")
        os.remove(path)


# ---------------------------------------------------------------------------
def get_hash(alg, file):
    assert alg in ("md5", "sha1", "sha256", "sha512")
    eprint(f"Computing {alg} hash of {str(file):<50}... ", end="")
    digest = from_cmd(alg+"sum", file).split(" ")[0]
    eprint(digest)
    return digest


# ---------------------------------------------------------------------------
g_lockfile = None  # Will delete in signal-handler, since @atexit handler
                   # does not fire when sync_files is interrupted when running
                   # in Docker container. GP-34570

def check_and_set_lockfile(lockfile_path) -> None:
    """ Check-or-create lockfile to verify that the directory is not currently in use by another process. """
    # This is similar to 'lockfile` unix command, but with some extra
    # functionality to detect and remove stale lockfiles.
    # We need it anyway since it's missing in busybox.

    hostname = os.uname()[1]
    wait_minutes = 0.125

    while lockfile_path.exists():
        (lf_pid, lf_hostname) = get_content(lockfile_path).strip().split(',')
        assert lf_hostname == hostname, f"The directory is locked by pid={lf_pid} on host {lf_hostname}.\nYou may need to remove the lockfile {lockfile_path} manually if that process is dead."

        if lf_pid in os.listdir('/proc') and os.getpgid(int(lf_pid)) != os.getpgid(os.getpid()):
            eprint(f"The directory is locked with {lockfile_path} by pid={lf_pid}. Waiting for {wait_minutes}m...")
            time.sleep(wait_minutes * 60)
            wait_minutes *= 2
        else:
            eprint(f"Clearing stale lockfile from pid {lf_pid}.", color=Color.BOLD)
            os.remove(lockfile_path)

    atexit.register(remove_if_exists, lockfile_path)

    global g_lockfile
    g_lockfile = lockfile_path

    try:
        os.makedirs(os.path.dirname(lockfile_path), exist_ok=True)
        with open(lockfile_path, 'x', encoding="ascii") as out:  # NB: opening in exclusive-creation mode
            print(os.getpid(), hostname, file=out, sep=',')
    except FileExistsError:
        raise
    except:
        # Other errors OK, e.g. don't have permissions
        pass


# ---------------------------------------------------------------------------
def make_manifest_file(dir, out_path):
    if dir[-1] == "/":
        dir = dir[:-1]

    assert dir[-1].isalnum(), dir

    files = [os.path.join(root, f) for root, _, files in os.walk(dir) for f in files]
    assert not any(out_path.endswith(file) for file in files), "Input directory contains output manifest path"

    for f in files:
        if any(s in f for s in ("mft", "manifest", "metalink")):
            eprint(f"Warning: the source directory contains a manifest file: {f}", color=Color.BOLD)

    data = {
        "version"    : 1,
        "totalFiles" : len(files),
        "timeStamp"  : datetime.now().isoformat(),
        "fileDetails":  [{
            "fileName"      : f[len(dir)+1:],      # Not a good name - actually sub-path in dir.
            "fileSize"      : os.path.getsize(f),  # NB: this field is string in current beta.
            "hashAlgorithm" : "md5",
            "hashValue"     : get_hash("md5", f)
        } for f in files ]
    }

    with open(out_path, "w", encoding="UTF-8") as out:
        json.dump(data, out, indent=4)


# ---------------------------------------------------------------------------
# NB: this is unused for now: maybe will switch from json-manifest to metalinks.
def make_metalink_file(dir, out_path):
    """
    Standand-based alternative to json-manifest.
    https://www.ietf.org/archive/id/draft-bryan-metalink-00.html
    """

    assert dir[-1].isalnum(), dir
    files = [os.path.join(root, f) for root, _, files in os.walk(dir) for f in files]
    assert not any(out_path.endswith(file) for file in files), "Input directory contains output metalink path"

    for f in files:
        if any(s in f for s in ("mft", "manifest", "metalink")):
            eprint(f"Warning: the source directory contains a manifest file: {f}", color=Color.BOLD)

    def make_file_elem(path, size, hash_digest):
        return f"""
        <file name="{path}">
            <size>{size}</size>
            <verification>
                <hash type="md5">{hash_digest}</hash>
            </verification>
            <resources>
                <url>{path}</url>
            </resources>
        </file>"""

    files_xml = "".join([make_file_elem(f[len(dir)+1:], os.path.getsize(f), get_hash("md5", f)) for f in files])

    xml = f"""<?xml version="1.0" encoding="UTF-8" ?>
<metalink version="3.0">
    <published>{datetime.now().isoformat()}</published>
    <files>{files_xml}
    </files>
</metalink>
"""

    with open(out_path, "w", encoding="ascii") as f:
        print(xml, file=f)


# ---------------------------------------------------------------------------
ManifestItem = namedtuple("ManifestItem", "path size hash_type hash_digest")
Manifest     = namedtuple("Manifest", "url text timestamp mtime files")


def get_mft(url) -> Manifest:
    """ Parse json-manifest or metalink """

    # If filepath-like, OK not to have a manifest
    if not url or (":/" not in str(url) and not Path(url).exists()):
        return None

    content = get_content(url)
    timestamp = None

    if content.startswith("<?xml"):  # metalink
        root  = ET.fromstring(content)
        assert root.tag == "metalink"
        timestamp = datetime.fromisoformat(root.findall("./published").text)
        files = [
            ManifestItem(
                f.find("resources/url").text,
                int(f.find("size").text),
                f.find("verification/hash").attrib["type"],
                f.find("verification/hash").text,
            ) for f in root.findall("./files/file")
        ]
    else:  # json-manifest
        data = json.loads(content)
        timestamp = datetime.fromisoformat(data["timeStamp"])
        files = [
            ManifestItem(f["fileName"], int(f["fileSize"]), f["hashAlgorithm"], f["hashValue"])
            for f in data["fileDetails"]
        ]

    assert all("://" not in f.path for f in files), "Expecting relative URLs only."

    # NB: mtime is the OS timestamp of the local manifest-file.
    # We'll put the manifest file in the output directory last after all files are
    # downloaded, and will compare file mtime against the manifest mtime to check whether
    # a file has been touched since it was downloaded.
    mtime = os.path.getmtime(url) if Path(url).exists() else None
    return Manifest(url, content, timestamp, mtime, files)


# ---------------------------------------------------------------------------
def file_integrity_ok(mi: ManifestItem, file_path, verify_hashes=False, verbose=False) -> bool:
    msg = (
             "file not found."    if not file_path.is_file()
        else "file-size changed." if os.path.getsize(file_path) != mi.size
        else "checksum changed."  if verify_hashes and mi.hash_digest != get_hash(mi.hash_type, file_path)
        else None
    )

    if verbose and msg:
        eprint(file_path, "-", msg, color=Color.RED)

    return not msg


# ---------------------------------------------------------------------------
def dir_integrity_ok(mft, dir, verify_hashes=False, verbose=False) -> bool:
    # NB: using sum() instead of all() to report all broken, not just first.
    return dir.is_dir() and len(mft.files) == sum(
        file_integrity_ok(mi, dir / mi.path, verify_hashes, verbose)
        for mi in mft.files
    )


# ---------------------------------------------------------------------------
def check_space(mft, dir):
    if not dir.exists():
        eprint(f"{dir} does not exist yet. Not checking disk-space.")
        return

    _, _, free     = shutil.disk_usage(dir)
    existing_paths = [dir / (f.path + sfx) for f in mft.files for sfx in ("", ".part")]
    existing_size  = sum((os.path.getsize(p) for p in existing_paths if p.exists()))
    incoming_size  = sum((f.size for f in mft.files))
    delta          = (incoming_size - existing_size)

    eprint("Space check: Available:", as_readable(free),
                       "; Existing:", as_readable(existing_size),
                       "; Incoming:", as_readable(incoming_size),
                          "; Delta:", as_readable(delta),
            sep = "",
    )
    assert free > delta, "Not enough free space available."


# ---------------------------------------------------------------------------
def download_file_with_aria(url, file_path):

    aria_config = file_path.with_suffix(".aria2_config")
    with open(aria_config, "w", encoding="ascii") as f:
        f.write("""
            file-allocation=none
            check-certificate=false
            allow-overwrite=true
            auto-file-renaming=false
            max-tries=5
            max-connection-per-server=5
            max-concurrent-downloads=5
            split=5
            console-log-level=warn
        """)

    subprocess.run(
        ["aria2c", f"--conf-path={aria_config}", f"--dir={file_path.parent}", f"--out={file_path.name}", url],
        check=True,
    )

    assert not file_path.with_suffix(".aria2").exists()
    assert file_path.exists()
    os.remove(aria_config)


# -----------------------------------------------------------------------
# Will compute hash on-the-fly, asynchronously in a singleton pool, as otherwise it becomes perf-bottleneck.
class async_hash:
    def __init__(self, hasher = hashlib.md5()):
        self.hash = hasher
        self.executor = concurrent.futures.ThreadPoolExecutor(1)
        self.fut = self.executor.submit(lambda: None)

    def update(self, chunk):
        self.fut.result()
        self.fut = self.executor.submit(self.hash.update, chunk)

    def hexdigest(self):
        self.fut.result()
        return self.hash.hexdigest()


# ---------------------------------------------------------------------------
# Context: rclone does not support resumable downloads natively,
# but rclone cat, fetching a chunk of file, is necessary and
# sufficient to implement parallelized resumable download of a large file.
def download_file_with_rclone(
    src_file,
    dst_file,
    chunk_size=100*1024*1024,
    num_streams=5
) -> str:  # md5 digest or None (can only compute on-the-fly if not resuming)
    """
    Copy file in chunks, in parallel.
    Resume partial transfers.
    Guard against source-file changes.
    """

    src_meta  = json.loads(from_cmd("rclone", "lsjson", src_file))[0]
    src_size  = src_meta["Size"]
    src_mtime = datetime.fromisoformat(src_meta["ModTime"][0:19])  # fromisoformat can't parse full ISO strings, so truncating
    dst_size  = 0 if not os.path.exists(dst_file) else os.path.getsize(dst_file)

    if os.path.exists(dst_file) and dst_size == src_size:
        return None # Download already complete - nothing to do

    # -----------------------------------------------------------------------
    def get_chunk(offset, this_chunk_size=None):  # -> (bytes, elapsed-time). Will be executed in parallel.
        if this_chunk_size is None:
            this_chunk_size = min(src_size, offset + chunk_size) - offset

        start_time = time.time()
        args = ["rclone", "cat", src_file, f"--offset={offset}", f"--count={this_chunk_size}"]
        # NB: check=False because rclone always errors-out when reading a chunk, at least on FTP
        buf = subprocess.run(args, check=False, capture_output=True).stdout
        assert this_chunk_size == len(buf), (this_chunk_size, len(buf))
        return (buf, time.time() - start_time)

    # TODO: add retry-wrapper over get_chunk?

    # -----------------------------------------------------------------------
    def get_bytes(file, offset, size):
        with open(file, "rb", encoding=None) as f:
            f.seek(offset)
            return f.read(size)

    # -----------------------------------------------------------------------
    def remove_dst_if_cant_resume():
        # Will check that tail of partially downloaded file
        # contains same bytes as in corresponding range in src.
        # We'll still check integrity downstream, but we want to bail early
        # if we can tell a-priori that resuming is pointless
        # (e.g. partial download is corrupted or src changed)
        nonlocal dst_size
        comp_offset = max(0, dst_size - 10*1024*1024)
        comp_size = dst_size - comp_offset
        can_resume = (
            0 < dst_size <= src_size   # otherwise src must have changed
            and src_mtime < datetime.fromtimestamp(os.path.getmtime(dst_file))
            and get_chunk(comp_offset, comp_size)[0] == get_bytes(dst_file, comp_offset, comp_size)
        )

        if can_resume:
            eprint(f"Resuming download of {src_file}.", color=Color.BOLD)
        elif os.path.exists(dst_file):
            eprint(f"Restarting download of {src_file} from the beginning.", color=Color.BOLD)
            remove_if_exists(dst_file)
            dst_size = 0

    remove_dst_if_cant_resume()

    # -----------------------------------------------------------------------
    start_time = time.time()
    starting_dst_size = dst_size
    md5 = async_hash(hashlib.md5())

    def write_chunk(chunk, out, elapsed):  # chunk of bytes, output file, and time it took to fetch the chunk
        md5.update(chunk)
        out.write(chunk)


        # stats and progress-bar
        avg_rate = int((dst_size - starting_dst_size)/(time.time() - start_time))
        prefix = "".join([
            Color.RESET,
            as_readable(dst_size), "/", as_readable(src_size),
            "(", str(int(dst_size*100/src_size)), "%) ",
            "CR:", as_readable(int(num_streams * chunk_size / elapsed)), "/s ",
            "AR:", as_readable(avg_rate), "/s ",
            "ETA:", str(timedelta(seconds=int((src_size - dst_size) / avg_rate))),
            " " * 8
        ])

        prog_bar = ""
        # prefix = prefix[:len(prefix) // 8 * 8]  # pad with spaces cap to mod-8 to deal with length-wobble.
        # max_prog_bar_width = (os.get_terminal_size().columns - len(prefix) + 6)
        # if max_prog_bar_width > 2:
        #    prog_bar_width = int(max_prog_bar_width * dst_size / src_size)
        #    pad_width = max_prog_bar_width - prog_bar_width - 1
        #    prog_bar = "[" + "="*prog_bar_width + " "*pad_width+ "]"

        eprint(prefix, prog_bar, end="\r")


    # -----------------------------------------------------------------------
    # map chunk-range -> chunk-bytes using a thread pool, and write them serially in this thread.
    # NB: can use multiprocessing.Pool.imap instead, but that's more overhead churning bytes.
    with concurrent.futures.ThreadPoolExecutor(num_streams) as pool:
        with open(dst_file, "ab", encoding=None) as out:  # in append-mode
            for (chunk, elapsed) in pool.map(get_chunk, range(dst_size, src_size, chunk_size)):
                dst_size += len(chunk)
                write_chunk(chunk, out, elapsed)
    eprint("")
    assert src_size == dst_size, (src_size, dst_size)

    return md5.hexdigest()


# ---------------------------------------------------------------------------
def transfer_file(mi: ManifestItem, src_mft_dir: str, out_dir: Path):
    """ Dispatch to appropriate downloader depending on url-type """

    file_path = out_dir / mi.path
    url = f"{src_mft_dir}/{mi.path}"

    # Make subdirs as necessary.
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # We'll be downloading to filename.part, and then replace the original (if exists),
    # but if it's "large" we delete the original to not run out of space.
    if file_path.exists() and os.path.getsize(file_path) > 1024**3:
        os.remove(file_path)

    is_fs = ":/" not in url
    tmp_file_path = file_path.parent / (file_path.name + ".part")
    md5 = None

    # -----------------------------------------------------------------------
    if is_fs:
        eprint(f"Copying {url} to {tmp_file_path}...")
        shutil.copy2(url, tmp_file_path)  # TODO: make resumable as with rclone-based implementation?

    elif url.startswith("gs://"):
        assert shutil.which("gsutil")
        subprocess.run(
            [
                "gsutil",
                "-o", "GSUtil:parallel_composite_upload_threshold=100M",
                "-o", "GSUtil:parallel_thread_count=2",
                "-o", "GSUtil:parallel_process_count=8",
                "-o", "GSUtil:sliced_object_download_max_components=8",
                "-m", "cp", "-L", str(tmp_file_path) + ".log",
                url, str(tmp_file_path),
            ],
            check=True
        )
        os.remove(str(tmp_file_path) + ".log")  # used by gsutil cp for resuming

    elif "://" in url and shutil.which("aria2c"):
        download_file_with_aria(url, tmp_file_path)

    elif "://" in url:  # fallback on curl if aria is not available
        subprocess.run(["curl", "-L", "-C", "-", "--retry", "5", "-o", tmp_file_path, url], check=True)

    else:
        assert ":/" in url  # my-rclone-remote:/path/to/file form
        assert shutil.which("rclone")
        md5 = download_file_with_rclone(url, tmp_file_path)

    # If computed hash on-the-fly, will not need to check it separately.
    hash_ok = is_fs or (md5 and md5 == mi.hash_digest)
    assert file_integrity_ok(mi, tmp_file_path, verify_hashes=(not hash_ok), verbose=True), "Integrity check failed."
    os.replace(tmp_file_path, file_path)

# ---------------------------------------------------------------------------
def can_rename(path) -> bool:
    # NB: For renaming dir to dir.in_progress during the download,
    # or for creating a lockfile in the parent dir, we can't do that
    # from within a container when dir is a mount-point.
    #
    # os.ismount() is not reliable and returns false in this scenario, eg.
    #
    # >> mkdir -p /tmp/foo/bar;
    # >> singularity exec --bind=/tmp/foo/bar fcs_img.sif python -c "import os; print(os.path.ismount('/tmp/foo/bar'));"
    # False
    #
    # So instead will check as follows
    try:
        os.rename(str(path), str(path) + "_tmp")
        os.rename(str(path) + "_tmp", str(path))
        return True
    except:
        return False


# ---------------------------------------------------------------------------
def main() -> None:

    # GP-34570
    # To propagate signals to child processes:
    # https://stackoverflow.com/questions/67823770/how-to-propagate-sigterm-to-children-created-via-subprocess
    def signal_handler(signo, _):
        remove_if_exists(g_lockfile)
        time.sleep(2) # give child processes (e.g. aria2c) a chance to wrap-up, e.g. to write a checkpoint for subsequent resume.
        raise SystemExit(f"\nSignal: {signo} ({signal.strsignal(signo)})")

    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGQUIT):
        signal.signal(sig, signal_handler)

    eprint("===============================================================================")

    args = parse_args()

    if args.mft.startswith("file://"):
        args.mft = args.mft[7:]

    if args.command == "make-manifest":
        make_manifest_file(args.dir, args.mft)
        sys.exit(0)

    src_mft       = get_mft(args.mft)
    src_mft_dir   = os.path.dirname(args.mft)
    src_mft_name  = os.path.basename(args.mft)
    dir           = Path(args.dir)
    assert dir not in (Path.home(), Path("."))

    is_dry_run    = args.command == "check" or dir in Path(args.mft).parents
    lockfile_path = Path(dir).with_suffix(".lockfile")

    if not is_dry_run:
        # NB: this needs to be done early, before accessing local_mft,
        # that will be deposited in dir by the process already running
        check_and_set_lockfile(lockfile_path)

    local_mft = get_mft(dir / src_mft_name)
    assert not dir.is_file(), f"{dir} is a file."

    if src_mft and local_mft and src_mft.timestamp != local_mft.timestamp:
        eprint(f"Remote manifest timestamp: {src_mft.timestamp}. Local manifest timestamp: {local_mft.timestamp}", color=Color.BOLD)

    if local_mft and local_mft.files == src_mft.files and dir_integrity_ok(src_mft, dir):
        eprint(f"{dir} is up-to-date with {src_mft_dir}.", color=Color.BOLD)
        sys.exit(0)

    assert args.command in ("get", "check")
    assert src_mft, "Could not access source-manifest."

    eprint(f"Source:      {src_mft_dir}")
    eprint(f"Destination: {dir}")

    if any((args.mft.startswith(s) for s in ["http://", "https://", "ftp://"])) and not shutil.which("aria2c"):
        eprint("Warning: aria2c is not accessible - will use curl instead (may be much slower).", color=Color.BOLD)

    # -----------------------------------------------------------------------
    # Prepare directory and check space.
    # To prevent leaving dir in invalid state, suffix it with .in_progress until done.
    # NB: it may already exist from a failed execution, in which case we'll reuse it.

    work_dir = dir.with_suffix(".in_progress") if can_rename(dir) else dir

    # -----------------------------------------------------------------------
    # Prepare work-dir.
    if work_dir.exists():
        eprint(f"Resuming failed transfer in {work_dir}...", color=Color.BOLD)
        assert local_mft is None
        local_mft = get_mft(work_dir / src_mft_name)

    elif is_dry_run:
        work_dir = dir
    elif dir.exists() and dir != work_dir:  # and is_in_place: # for now always-in-place
        os.rename(dir, work_dir)
    else:
        os.makedirs(work_dir)

    check_space(src_mft, work_dir)

    # -----------------------------------------------------------------------
    # Transfer each file in the manifest, smaller first
    src_mft.files.sort(key=lambda mi: mi.size)
    for mi in src_mft.files:
        eprint()
        file_path    = work_dir / mi.path
        in_local_mft = local_mft and (mi in local_mft.files)


        assert "../" not in mi.path
        assert work_dir in file_path.parents


        file_done = file_integrity_ok(
            mi,
            file_path,
            verbose=in_local_mft,  # if we have it and integrity-check fails - report.
            verify_hashes=not (in_local_mft and file_path.exists() and os.path.getmtime(file_path) < local_mft.mtime)
        )


        eprint("Skipping existing" if file_done else "Requires transfer:", f"{as_readable(mi.size)} {mi.path}", color=Color.RESET)
        if not file_done and not is_dry_run:
            transfer_file(mi, src_mft_dir, work_dir)

    # -----------------------------------------------------------------------
    # Final steps: put the manifest in the output dir and rename it back.
    if is_dry_run:
        return

    assert dir_integrity_ok(src_mft, work_dir)

    with open(work_dir / src_mft_name, "w", encoding="ascii") as f:
        f.write(src_mft.text)

    if dir.exists() and dir != work_dir:
        shutil.rmtree(dir)

    remove_if_exists(lockfile_path)  # before renaming dir

    if dir != work_dir:
        os.rename(work_dir, dir)

    eprint("Done.", color=Color.BOLD)


if __name__ == "__main__":
    main()



#########################################################################################################
# Short tutorial:
#########################################################################################################

# disable warnings about "pointless" string statements below
# pylint: disable=W0105

"""
export MFT_URL=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest
# If URL is an ordinary URL, will use aria2c if available, otherwise curl.
# Otherwise, if URL is in the form ftp-gcs-rclone-remote:/genomes/TOOLS/FCS/database/test-only/test-only.manifest, will use rclone.
# Otherwise, will treat URL as local file path.


# Download db directly to shm:
./scripts/sync_files.py get --mft=$MFT_URL --dir=/dev/shm/gxdb/test-only

# Or, download db to disk, and then from disk to shm.
./scripts/sync_files.py get --mft=$MFT_URL --dir=./tmp/test-only
./scripts/sync_files.py get --mft=./tmp/test-only/test-only.manifest --dir=/dev/shm/gxdb/test-only

# To check if data changed upstream:
./scripts/sync_files.py check --mft=$MFT_URL --dir=/dev/shm/gxdb/test-only

# To update, run the commands again (they do nothing if files are in sync).
"""


#########################################################################################################
# Long(er) tutorial:
#########################################################################################################
"""
export DOWNLOAD_CMD="./scripts/sync_files.py --mft=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest --dir=tmp/test-only"
rm -rf tmp/test-only


$DOWNLOAD_CMD check             # See what's going to be transfered.
timeout 5s $DOWNLOAD_CMD get    # Start the download and kill it after 5 seconds.

ls tmp/test-only.in_progress/   # Observe that tmp/test-only has been renamed to tmp/test-only.in_progress during the download,
                                # such that at no time tmp/test-only is available in an invalid state.

$DOWNLOAD_CMD get               # Resume the download. Observe the message about resuming a failed transfer.
$DOWNLOAD_CMD get               # Observe the message that there's nothing to do.


# Corrupt a few files:
rm tmp/test-only/test-only.blast_div.tsv.gz                     # remove the file
perl -pi -e 's/build/BUILD/' tmp/test-only/test-only.meta.jsonl # same size, but content changed
echo foo > tmp/test-only/test-only.gxi                          # size changed
touch tmp/test-only/test-only.gxs                               # timestamp changed, but content didn't


# Check integrity using the local manifest.
# Observe the messages about missing/corrupted files.
./scripts/sync_files.py check --mft=tmp/test-only/test-only.manifest --dir=tmp/test-only

# Launch command in background; launch duplicate command in foreground.
# Observe messages from the duplicate command that the directory is locked

# and it's waiting for the other one to complete.
# After the background command completes the work and exits,
# the second command reports that everything is up-to-date and exits.
$DOWNLOAD_CMD get &
sleep 1
$DOWNLOAD_CMD get
wait
"""
