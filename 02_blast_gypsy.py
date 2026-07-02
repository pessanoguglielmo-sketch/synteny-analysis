"""
01_blast_gypsy.py

Filters Gypsy-family transposable elements from a TE annotation FASTA,
builds a BLAST nucleotide database from a target genome, and BLASTs the
Gypsy sequences against that genome.

Requirements:
- Biopython (pip install biopython)
- NCBI BLAST+ (makeblastdb, blastn) installed and available on PATH,
  or set the BLAST_BIN environment variable to the folder containing them.

Inputs (expected in the working directory):
- <genome_fa>   : target genome FASTA (e.g. DS_hapA.fasta)
- TE_GD.fa      : transposable element annotation FASTA containing Gypsy entries

Outputs:
- gypsy_only.fa           : filtered Gypsy-only FASTA
- genome_db4.*             : BLAST database files
- gypsy_vs_genomeDS.tsv    : BLAST hits (outfmt 6)
"""

from Bio import SeqIO
import subprocess
import os

# ---------- 1) Config ----------
genome_fa = "DS_hapA.fasta"  # target genome to BLAST against
te_annotation_fa = "TE_GD.fa"
blast_db_name = "genome_db4"
blast_output_tsv = "gypsy_vs_genomeDS.tsv"

# BLAST+ binaries: uses BLAST_BIN env var if set, otherwise assumes
# makeblastdb/blastn are on PATH.
blast_bin = os.environ.get("BLAST_BIN", "")
makeblastdb = os.path.join(blast_bin, "makeblastdb") if blast_bin else "makeblastdb"
blastn = os.path.join(blast_bin, "blastn") if blast_bin else "blastn"

# ---------- 2) Filter Gypsy elements ----------
gypsy_records = [
    rec for rec in SeqIO.parse(te_annotation_fa, "fasta")
    if "gypsy" in rec.description.lower()
]

SeqIO.write(gypsy_records, "gypsy_only.fa", "fasta")
print(f"Wrote {len(gypsy_records)} Gypsy records to gypsy_only.fa")

# ---------- 3) Build BLAST database from target genome ----------
subprocess.run([
    makeblastdb,
    "-in", genome_fa,
    "-dbtype", "nucl",
    "-out", blast_db_name
], check=True)
print(f"BLAST database created: {blast_db_name}.*")

# ---------- 4) Run BLAST: Gypsy vs genome ----------
subprocess.run([
    blastn,
    "-query", "gypsy_only.fa",
    "-db", blast_db_name,
    "-evalue", "1e-5",
    "-outfmt", "6",
    "-out", blast_output_tsv
], check=True)
print(f"BLAST finished. Results in {blast_output_tsv}")
