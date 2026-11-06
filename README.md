# copyratio

Tell plasmid carriage apart from library carry-over, using read depth rather
than sequence comparison.

## The problem

A plasmid sequence can turn up in an isolate's assembly without that isolate
carrying it. Contaminating DNA from another sample in the same sequencing run
assembles perfectly well, and the resulting contig looks like any other.

Sequence-based methods cannot settle this, because they all ask the same
question: *is this sequence present?* Contamination and carriage both answer
yes. BLAST coverage, plasmid classifiers and assembly graph structure are all
variations on that one question, so they tend to agree with each other whether
or not they are right.

`copyratio` asks a different question: *how much of it is there, relative to
the cell it is supposed to be inside?*

## Why the ratio works

Inside a cell that genuinely carries a plasmid, the plasmid is present at
roughly its copy number — commonly several copies per chromosome. So if you
sequence that culture, the plasmid comes out **at or above** chromosome depth.

If instead a trace of foreign DNA entered the library, the plasmid arrives at
whatever fraction that contamination represents. That is typically a few
percent, so the plasmid sits **far below** the host's own chromosome.

The ratio is therefore an estimate of copies per cell, and it separates the two
cases cleanly even when coverage breadth does not.

**The denominator has to be the isolate's own assembly.** Not a shared
reference, not a related strain. Genome sizes differ several-fold across
species, and sequencing yield differs between libraries; both cancel out only
when the same reads are mapped to the same cell's own genome.

## What you need

- `bwa`
- `samtools`
- Python 3.7+
- For each isolate: its reads, and an assembly of that isolate
- The element you want to test, as a FASTA

## Install

```bash
git clone https://github.com/USERNAME/copyratio.git
cd copyratio
chmod +x copyratio.py
```

Or with conda for the dependencies:

```bash
conda create -n copyratio -c bioconda bwa samtools python=3.11
conda activate copyratio
```

## Use

One isolate:

```bash
./copyratio.py \
    --element plasmid.fasta \
    --host ASMA-1664_assembly.fna \
    --r1 ASMA-1664_R1.fastq.gz \
    --r2 ASMA-1664_R2.fastq.gz \
    --sample ASMA-1664 \
    -o result.tsv
```

Many isolates — write a manifest with columns `sample`, `r1`, `r2`, `host`:

```
sample      r1                  r2                  host
ASMA-1664   1664_R1.fq.gz       1664_R2.fq.gz       1664_assembly.fna
ASMA-1632   1632_R1.fq.gz       1632_R2.fq.gz       1632_assembly.fna
ASMA-1667   1667_R1.fq.gz       1667_R2.fq.gz       1667_assembly.fna
```

```bash
./copyratio.py --element plasmid.fasta --manifest samples.tsv -o results.tsv -t 8
```

## Output

| column | meaning |
| --- | --- |
| `element_depth` | mean depth across the element, zero positions included |
| `element_breadth` | percent of the element with any coverage |
| `host_depth` | mean depth across that isolate's own assembly |
| `ratio` | `element_depth / host_depth` — roughly copies per cell |
| `fold_below` | `1 / ratio`, for readability when the ratio is small |
| `call` | `carriage`, `carry-over`, `ambiguous` or `undetermined` |

Calls come from `--carriage-ratio` (default 0.5) and `--carryover-ratio`
(default 0.1). Anything between them is `ambiguous` and worth looking at by
hand. `undetermined` means host depth fell below `--min-host-depth`, where the
ratio stops being stable.

**The thresholds are conventions, not constants.** Read the ratio, not the
call.

## A worked example

Three non-*Pseudomonas* isolates from a cystic fibrosis cohort appeared to
carry a *P. aeruginosa* resistance plasmid, at 54–95% BLAST coverage. One of
them came from a different patient, which would have made it cross-species and
cross-patient transfer.

| isolate | species | element depth | host depth | ratio |
| --- | --- | --- | --- | --- |
| ASMA-1664 | *Streptococcus salivarius* | 5.1× | 160× | 0.032 |
| ASMA-1632 | *Streptococcus* sp. | 3.6× | 244× | 0.015 |
| ASMA-1662 | *Bacillus cereus* | 1.6× | 86× | 0.019 |
| ASMA-1667 | *P. aeruginosa* (true carrier) | 286× | 61× | **4.7** |

The three candidates sit 31–68× below their own chromosomes. The genuine
carrier sits 4.7× above its own, which is what a multi-copy plasmid should do.
All three were library carry-over.

Mapping the same reads against the *P. aeruginosa* chromosome afterwards
confirmed it: the chromosome was present too, at 0.6–1.6×, spread over 22–67%
of its 7 Mb. The contaminating DNA was a whole genome at its natural plasmid
copy number, not a loose plasmid.

## Limits

**It needs a decent assembly of the host.** The denominator is only meaningful
if the assembly represents that isolate's genome. A heavily contaminated
assembly contains the contaminant too, which drags the denominator in the wrong
direction. Check your host assemblies first.

**It is unreliable at low depth.** Below roughly 10× on the host, the ratio
swings on a handful of reads. `--min-host-depth` enforces this.

**Single-copy and integrated elements will sit near 1.0, not above it.** The
0.5 default assumes nothing about copy number beyond "comparable to the
chromosome." If you are testing an element you expect to be chromosomally
integrated, carriage looks like ~1.0, and you should set the thresholds
accordingly.

**It says whether a cell carries the element. It does not say how
contamination happened.** Index hopping, library cross-contamination and
plate-level carry-over all produce the same signature here. Distinguishing them
needs sequencing batch metadata, which this tool does not look at.

**Repetitive elements can inflate the numerator.** If the element shares
sequence with something genuinely in the host genome — a transposase family,
say — those reads map to the element and raise its depth. Mapping quality
filtering (`--mapq`, default 20) helps but does not eliminate this. If the
element's coverage is concentrated in one region rather than spread across it,
treat the ratio with caution and look at the per-position profile.

## Licence

MIT
