#!/usr/bin/env python3
"""
copyratio — tell plasmid carriage apart from library carry-over.

A sequence can be present in an isolate's reads without the cell carrying it.
This computes how deeply the element is sequenced relative to that same
isolate's own genome. Carriage puts the two at comparable depth; carry-over
leaves the element far below.

See README.md for why the host's own assembly has to be the denominator.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile


def need(tool):
    if shutil.which(tool) is None:
        sys.exit("error: %s not found on PATH" % tool)


def index_if_absent(fasta, workdir):
    """bwa needs an index next to the fasta. If we can't write there, copy it."""
    if all(os.path.exists(fasta + ext) for ext in (".bwt", ".pac", ".sa")):
        return fasta
    target = fasta
    if not os.access(os.path.dirname(os.path.abspath(fasta)) or ".", os.W_OK):
        target = os.path.join(workdir, os.path.basename(fasta))
        shutil.copyfile(fasta, target)
    subprocess.run(["bwa", "index", target], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return target


def depth_stats(reference, r1, r2, threads, mapq, workdir, keep_bam=None):
    """Map reads to a reference and return (mean depth, breadth %, length)."""
    ref = index_if_absent(reference, workdir)
    bam = keep_bam or os.path.join(workdir, "aln.bam")

    bwa = subprocess.Popen(
        ["bwa", "mem", "-t", str(threads), ref, r1] + ([r2] if r2 else []),
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    sort = subprocess.Popen(
        ["samtools", "sort", "-@", str(max(1, threads // 2)), "-o", bam],
        stdin=bwa.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    bwa.stdout.close()
    sort.communicate()
    if sort.returncode != 0:
        raise RuntimeError("alignment failed for %s" % reference)
    subprocess.run(["samtools", "index", bam], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # -a emits zero-coverage positions; without it the mean is inflated.
    out = subprocess.run(["samtools", "depth", "-a", "-Q", str(mapq), bam],
                         check=True, stdout=subprocess.PIPE,
                         text=True).stdout

    total = covered = positions = 0
    for line in out.splitlines():
        d = int(line.rsplit("\t", 1)[1])
        total += d
        positions += 1
        if d > 0:
            covered += 1
    if positions == 0:
        return 0.0, 0.0, 0
    return total / positions, covered / positions * 100, positions


def run_one(sample, r1, r2, element, host, args, out_writer):
    with tempfile.TemporaryDirectory(dir=args.tmpdir) as work:
        e_depth, e_breadth, e_len = depth_stats(
            element, r1, r2, args.threads, args.mapq, work)
        h_depth, h_breadth, h_len = depth_stats(
            host, r1, r2, args.threads, args.mapq, work)

    ratio = (e_depth / h_depth) if h_depth > 0 else float("nan")
    if h_depth < args.min_host_depth:
        call = "undetermined"
        note = "host depth %.1fx below --min-host-depth" % h_depth
    elif ratio >= args.carriage_ratio:
        call = "carriage"
        note = ""
    elif ratio <= args.carryover_ratio:
        call = "carry-over"
        note = ""
    else:
        call = "ambiguous"
        note = "ratio falls between the two thresholds"

    row = dict(
        sample=sample,
        element_depth="%.2f" % e_depth,
        element_breadth="%.1f" % e_breadth,
        element_bp=e_len,
        host_depth="%.2f" % h_depth,
        host_breadth="%.1f" % h_breadth,
        host_bp=h_len,
        ratio="%.3f" % ratio,
        fold_below=("%.0f" % (1 / ratio)) if ratio > 0 else "",
        call=call,
        note=note,
    )
    out_writer.writerow(row)
    sys.stderr.write(
        "  %-18s element %7.2fx (%.0f%% covered)   host %7.2fx   "
        "ratio %.3f   %s\n"
        % (sample, e_depth, e_breadth, h_depth, ratio, call))
    return row


FIELDS = ["sample", "element_depth", "element_breadth", "element_bp",
          "host_depth", "host_breadth", "host_bp", "ratio", "fold_below",
          "call", "note"]


def main():
    p = argparse.ArgumentParser(
        description="Distinguish plasmid carriage from library carry-over "
                    "by read-depth ratio.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
single sample:
  copyratio.py --element plasmid.fa --host asm.fa \\
      --r1 s_R1.fq.gz --r2 s_R2.fq.gz --sample S1 -o out.tsv

batch (manifest is a TSV with columns sample, r1, r2, host):
  copyratio.py --element plasmid.fa --manifest samples.tsv -o out.tsv
""")
    p.add_argument("--element", required=True,
                   help="FASTA of the sequence being tested (plasmid, "
                        "integron, any element)")
    p.add_argument("--manifest",
                   help="TSV with columns: sample, r1, r2, host")
    p.add_argument("--sample", help="sample name (single-sample mode)")
    p.add_argument("--r1", help="forward reads (single-sample mode)")
    p.add_argument("--r2", help="reverse reads; omit for single-end")
    p.add_argument("--host", help="this isolate's own assembly "
                                  "(single-sample mode)")
    p.add_argument("-o", "--out", default="-",
                   help="output TSV, '-' for stdout (default)")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--mapq", type=int, default=20,
                   help="minimum mapping quality (default 20)")
    p.add_argument("--carriage-ratio", type=float, default=0.5,
                   help="ratio at or above which to call carriage "
                        "(default 0.5)")
    p.add_argument("--carryover-ratio", type=float, default=0.1,
                   help="ratio at or below which to call carry-over "
                        "(default 0.1)")
    p.add_argument("--min-host-depth", type=float, default=10.0,
                   help="below this the ratio is unstable (default 10)")
    p.add_argument("--tmpdir", default=None)
    args = p.parse_args()

    need("bwa")
    need("samtools")

    jobs = []
    if args.manifest:
        with open(args.manifest) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                missing = [c for c in ("sample", "r1", "host") if not rec.get(c)]
                if missing:
                    sys.exit("manifest row missing %s: %r" % (missing, rec))
                jobs.append((rec["sample"], rec["r1"],
                             rec.get("r2") or None, rec["host"]))
    else:
        if not (args.sample and args.r1 and args.host):
            p.error("need --manifest, or --sample/--r1/--host together")
        jobs.append((args.sample, args.r1, args.r2, args.host))

    for s, r1, r2, host in jobs:
        for f in [r1, host] + ([r2] if r2 else []):
            if not os.path.exists(f):
                sys.exit("error: no such file: %s" % f)

    fh = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()

    sys.stderr.write("element: %s\n%d sample(s)\n\n" % (args.element, len(jobs)))
    for s, r1, r2, host in jobs:
        try:
            run_one(s, r1, r2, args.element, host, args, w)
        except Exception as exc:
            sys.stderr.write("  %-18s FAILED: %s\n" % (s, exc))
            w.writerow(dict(sample=s, call="failed", note=str(exc)))
        fh.flush()

    if fh is not sys.stdout:
        fh.close()
        sys.stderr.write("\nwrote %s\n" % args.out)


if __name__ == "__main__":
    main()
