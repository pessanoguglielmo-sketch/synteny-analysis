#!/bin/bash
# 01a_synteny_alignment.slurm
#
# SLURM job for the KU Leuven VSC HPC (Genius cluster).
# Aligns a reference genome (GD haplome 1) against several other Maleae
# genomes/haplomes using minimap2, sorts/indexes the resulting BAMs,
# then runs SyRI to detect structural rearrangements and plotsr to
# visualise them.
#
# Requirements:
# - minimap2, samtools (module system or conda)
# - conda envs: one containing samtools/seqtk (here "blastn"),
#   one containing syri + plotsr (here "syri_env")
# - Genome FASTA files, pre-renamed so chromosome IDs match across genomes
#   (see chrXX naming convention used elsewhere in this project)
#
# Usage: adjust the SBATCH resource requests, --account, --mail-user,
# working directory (WORKDIR) and genome paths below for your own
# cluster allocation before submitting with `sbatch`.

#SBATCH --job-name=synteny_all
#SBATCH --account=<your_project_account>
#SBATCH --clusters=genius
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --time=03:00:00
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH -o %j_out_synteny_all.out
#SBATCH -e %j_err_synteny_all.out

date

# ---------- Activate environment (samtools, seqtk) ----------
source ~/.bashrc
source activate blastn

# ---------- Load modules ----------
module purge
module load cluster/genius/bigmem
module load minimap2/2.26-GCCcore-12.3.0

# ---------- Working directory ----------
WORKDIR=./alignments
cd "$WORKDIR"

# ---------- Genome paths (relative to WORKDIR's parent) ----------
genomeA=../GD_haplome1_renamed.fasta
genomeB=../GD_haplome2_renamed.fasta
genomeC=../dA_haplome1_renamed.fasta
genomeD=../dA_haplome2_renamed.fasta
genomeE=../DS_hapA_renamed.fasta
genomeF=../DS_hapB_renamed.fasta
genomeG=../MF_hap1_renamed.fasta
genomeH=../MF_hap2_renamed.fasta
genomeI=../loquat_genome_renamed.fasta
genomeJ=../crataegus_genome_chr.fasta

######### Align the genomes ##################

# apple haplomes
#minimap2 -x asm5 -t 10 --eqx $genomeA $genomeB > A_B.paf
#minimap2 -ax asm5 -t 10 --eqx $genomeA $genomeB > A_B.sam
#samtools view -bS A_B.sam | samtools sort -o A_B_sort.bam -
#samtools index A_B_sort.bam

# pyrus haplome1 to apple haplome1
#minimap2 -x asm5 -t 10 --eqx $genomeA $genomeC > A_C.paf
#minimap2 -ax asm5 -t 10 --eqx $genomeA $genomeC > A_C.sam
#samtools view -bS A_C.sam | samtools sort -o A_C_sort.bam -
#samtools index A_C_sort.bam

# pyrus haplome2 to apple haplome1
#minimap2 -x asm5 -t 10 --eqx $genomeA $genomeD > A_D.paf
#minimap2 -ax asm5 -t 10 --eqx $genomeA $genomeD > A_D.sam
#samtools view -bS A_D.sam | samtools sort -o A_D_sort.bam -
#samtools index A_D_sort.bam

# crataegus genome to apple haplome1
minimap2 -x asm5 -t 10 --eqx $genomeA $genomeJ > A_J.paf

# Pyrus pyrifolia haplomeA to apple haplome1
#minimap2 -x asm5 -t 10 --eqx $genomeA $genomeE > A_E.paf
#minimap2 -ax asm5 -t 10 --eqx $genomeA $genomeE > A_E.sam
#samtools view -bS A_E.sam | samtools sort -o A_E_sort.bam -
#samtools index A_E_sort.bam

# Pyrus pyrifolia haplomeB to apple haplome1
#minimap2 -x asm5 -t 10 --eqx $genomeA $genomeF > A_F.paf
#minimap2 -ax asm5 -t 10 --eqx $genomeA $genomeF > A_F.sam
#samtools view -bS A_F.sam | samtools sort -o A_F_sort.bam -
#samtools index A_F_sort.bam

# Malus fusca haplome 1 to apple haplome1
minimap2 -x asm5 -t 10 --eqx $genomeA $genomeG > A_G.paf
minimap2 -ax asm5 -t 10 --eqx $genomeA $genomeG > A_G.sam
samtools view -bS A_G.sam | samtools sort -o A_G_sort.bam -
samtools index A_G_sort.bam

# Malus fusca haplome 2 to apple haplome1
minimap2 -x asm5 -t 10 --eqx $genomeA $genomeH > A_H.paf
minimap2 -ax asm5 -t 10 --eqx $genomeA $genomeH > A_H.sam
samtools view -bS A_H.sam | samtools sort -o A_H_sort.bam -
samtools index A_H_sort.bam

# loquat genome to apple haplome1
minimap2 -x asm5 -t 10 --eqx $genomeA $genomeI > A_I.paf


######### Perform synteny analysis and plot ##################

# ---------- Activate environment (syri, plotsr) ----------
source ~/.bashrc
source activate syri_env

# run syri
#syri -c A_B_sort.bam -r $genomeA -q $genomeB -F B --prefix A_B
#syri -c A_C_sort.bam -r $genomeA -q $genomeC -F B --prefix A_C
#syri -c A_D_sort.bam -r $genomeA -q $genomeD -F B --prefix A_D
#syri -c A_E_sort.bam -r $genomeA -q $genomeE -F B --prefix A_E
#syri -c A_F_sort.bam -r $genomeA -q $genomeF -F B --prefix A_F
syri -c A_G_sort.bam -r $genomeA -q $genomeG -F B --prefix A_G
syri -c A_H_sort.bam -r $genomeA -q $genomeH -F B --prefix A_H

# run plotsr
#plotsr --sr A_Bsyri.out --genomes ../genomes_AB.txt -o plot_AB.png
#plotsr --sr A_Csyri.out --genomes ../genomes_AC.txt -o plot_AC.png
#plotsr --sr A_Dsyri.out --genomes ../genomes_AD.txt -o plot_AD.png
plotsr --sr A_Esyri.out --genomes ../genomes_AE.txt -o plot_AE.png
plotsr --sr A_Fsyri.out --genomes ../genomes_AF.txt -o plot_AF.png
plotsr --sr A_Gsyri.out --genomes ../genomes_AG.txt -o plot_AG.png
plotsr --sr A_Hsyri.out --genomes ../genomes_AH.txt -o plot_AH.png
