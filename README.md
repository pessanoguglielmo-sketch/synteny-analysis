# Apple–Pear Synteny Analysis

Comparative genomics pipeline analyzing genome structure and synteny between
apple (*Malus*) and related Maleae genomes (pear, *Pyrus*; hawthorn,
*Crataegus*; loquat, *Eriobotrya*; *Malus fusca*), developed as part of an
MSc thesis on apple–pear hybridisation and genome structure at the Plant
Genomics & Crop Improvement Lab, KU Leuven.

The pipeline combines whole-genome alignment, Gypsy transposable element
(TE) density, and optimal chromosome-pair matching to characterize
large-scale structural conservation and rearrangement between genome pairs,
and to test whether Gypsy elements co-localize with synteny-disruption
sites.

## Pipeline overview

```
                ┌─────────────────────────────┐
                │   01a_synteny_alignment.slurm │
                │   minimap2 whole-genome        │
                │   alignment (PAF) + SyRI/plotsr│
                └───────────────┬─────────────┘
                                │
┌─────────────────────────────┐│
│  01b_blast_gypsy.py           ││
│  Filter Gypsy TEs, BLAST       ││
│  against target genome (TSV)   ││
└───────────────┬─────────────┘│
                │                │
                └───────┬────────┘
                        ▼
              ┌───────────────────────┐
              │  02_synteny_plot.R      │
              │  - Chromosome sanity    │
              │    checks                │
              │  - Optimal chromosome    │
              │    pairing (Hungarian    │
              │    algorithm)             │
              │  - Reciprocal synteny %   │
              │  - Per-chromosome-pair    │
              │    synteny + Gypsy-density│
              │    plots (gggenomes)      │
              └───────────────────────┘
```

`01a` and `01b` are independent and can run in either order (or in
parallel) — `01a` produces the whole-genome alignment (PAF), `01b`
produces the Gypsy-vs-genome BLAST hits. `02_synteny_plot.R` requires
outputs from both.

## Repository contents

| File | Description |
|---|---|
| `01a_synteny_alignment.slurm` | SLURM job (KU Leuven VSC Genius cluster) that aligns a reference genome against several Maleae genomes/haplomes with `minimap2`, sorts/indexes the resulting BAMs, then runs `SyRI` to call structural rearrangements and `plotsr` to visualise them. |
| `01b_blast_gypsy.py` | Filters Gypsy-family TE sequences from a transposable element annotation FASTA, builds a BLAST nucleotide database from the target genome, and BLASTs the Gypsy sequences against it. |
| `02_synteny_plot.R` | Parses Gypsy TE density (from GFF and BLAST TSV), reads genome FASTAs and the PAF alignment, finds the optimal one-to-one chromosome pairing that maximises synteny, computes strict reciprocal synteny, and generates per-chromosome synteny + Gypsy-density plots. |

## Methodology notes

**Why optimal chromosome pairing instead of naive chr01↔chr01 matching?**
Genomes that have undergone structural rearrangement (fusions, fissions,
translocations) don't necessarily have chromosome *N* in genome A
corresponding to chromosome *N* in genome B. `02_synteny_plot.R` builds a
full pairwise synteny-coverage matrix across all chromosome combinations,
then uses the Hungarian algorithm (`clue::solve_LSAP`) to find the
one-to-one assignment that maximises total synteny across all pairs. This
is a more defensible basis for downstream analysis (e.g. testing Gypsy
co-localization with synteny breakpoints) than assuming numeric
correspondence.

**Why "strict" reciprocal synteny?**
Overall synteny percentage is computed using *only* the chromosome pairs
selected by the optimal assignment above — not all pairwise alignment
hits. This avoids inflating the synteny estimate with spurious
cross-chromosome alignments that aren't part of the true one-to-one
chromosome correspondence.

**Chromosome ID sanity checks.**
Chromosome numbers are extracted from sequence IDs via regex
(`extract_chr_num()`), which can fail silently if header formats vary
unexpectedly. `check_chr_ids()` runs after every extraction step and:
- **stops** the script if two different sequence IDs map to the same
  chromosome number (a sign the regex is mis-parsing headers), reporting
  which IDs collided;
- **warns** if the number of distinct chromosomes found doesn't match the
  expected `N_CHRS`, without halting execution.

## Requirements

**Cluster / alignment (`01a`):**
- SLURM job scheduler (KU Leuven VSC Genius, or adapt for your own cluster)
- `minimap2`, `samtools`
- `syri`, `plotsr` (conda environment)

**BLAST (`01b`):**
- Python 3, [Biopython](https://biopython.org/)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`makeblastdb`, `blastn`)

**Plotting / synteny analysis (`02`):**
- R (≥ 4.0 recommended)
- CRAN packages: `ggplot2`, `dplyr`, `data.table`, `clue`, `scales`
- Bioconductor packages: `GenomicRanges`, `IRanges`
  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("GenomicRanges", "IRanges"))
  ```
- [`gggenomes`](https://github.com/thackl/gggenomes)

## Usage

1. **Prepare genomes**: rename chromosomes to a consistent format (e.g.
   `chr01`, `chr02`, ...) across all genomes being compared, so
   `extract_chr_num()` can parse them correctly.

2. **Run the alignment step** (`01a_synteny_alignment.slurm`): edit the
   `--account`, `--mail-user` (if desired), working directory, and genome
   paths for your own cluster allocation, then submit with:
   ```bash
   sbatch 01a_synteny_alignment.slurm
   ```

3. **Run the Gypsy BLAST step** (`01b_blast_gypsy.py`): edit the
   `genome_fa`, TE annotation filename, and output filenames at the top of
   the script for your genome pair, then run:
   ```bash
   python 01b_blast_gypsy.py
   ```
   (Set the `BLAST_BIN` environment variable if `makeblastdb`/`blastn`
   aren't on your `PATH`.)

4. **Run the synteny analysis and plotting step** (`02_synteny_plot.R`):
   edit the `USER INPUTS` section at the top (file paths, number of
   chromosomes, window size, alignment/BLAST filtering thresholds) to
   match your genome pair and the outputs from steps 1–2, then run:
   ```r
   source("02_synteny_plot.R")
   ```

## Outputs

Written to the `OUTDIR` folder configured in `02_synteny_plot.R`:

| File | Description |
|---|---|
| `all_pairwise_chromosome_synteny.csv` | Synteny coverage stats for every chromosome-pair combination considered. |
| `assigned_best_pairs.csv` | The optimal one-to-one chromosome assignment (Hungarian algorithm) with per-pair synteny stats. |
| `overall_reciprocal_synteny_strict.csv` | Genome-wide reciprocal synteny percentage, computed using only the assigned best chromosome pairs. |
| `synteny_<chrA>_vs_<chrB>.png` | One synteny + Gypsy-density plot per assigned chromosome pair. |

## Notes

- Genome FASTA and alignment/annotation files (FASTA, GFF, PAF, BAM, BLAST
  TSVs) are not included in this repository due to size; scripts assume
  these are provided locally.
- This is analysis code developed for a specific thesis project rather
  than a general-purpose packaged tool — file paths, genome labels, and
  thresholds in the `USER INPUTS` sections are meant to be edited per run
  rather than passed as command-line arguments.
