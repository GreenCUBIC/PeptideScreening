#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --gpus-per-node=4
#SBATCH --nodes=1
#SBATCH -p compute_full_node
#SBATCH --account=soscip-3-091
#SBATCH --output=/home/j/jrgreen/fcharih/scratch/peptide-screening/logs/predict_dscript-%x-%j.out

EMBEDDINGS=$HOME/scratch/peptide-screening/data/outputs/dscript/human_proteome_with_peptides.embed
PAIRS=$HOME/scratch/peptide-screening/data/inputs/pairs_to_score_pipr_without_header.tsv
OUTPUT_FILE=$HOME/scratch/peptide-screening/data/outputs/dscript/predictions.txt
MODEL=$HOME/scratch/peptide-screening/data/inputs/human_v1.sav

# Load the modules
module load anaconda3 gcc

# Load the virtual env
cd $HOME
source activate dscriptenv

dscript predict -o $OUTPUT_FILE -d 0 --model $MODEL --pairs $PAIRS --embeddings $EMBEDDINGS 
