#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --gpus-per-node=4
#SBATCH --nodes=1
#SBATCH -p compute_full_node
#SBATCH --account=soscip-3-091
#SBATCH --output=/home/j/jrgreen/fcharih/scratch/peptide-screening/logs/embed_dscript-%x-%j.out

INFILE=$HOME/scratch/peptide-screening/data/inputs/human_proteome_with_peptides.fasta
OUTFILE=$HOME/scratch/peptide-screening/data/outputs/human_proteome_with_peptides.embed

# Load the modules
module load anaconda3 gcc

# Load the virtual env
cd $HOME
source activate dscriptenv

dscript embed --seqs $INFILE --outfile $OUTFILE
