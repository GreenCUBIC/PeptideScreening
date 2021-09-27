#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --gpus-per-node=4
#SBATCH --nodes=1
#SBATCH -p compute_full_node
#SBATCH --account=soscip-3-091
#SBATCH --output=/home/j/jrgreen/fcharih/scratch/peptide-screening/logs/pipr-prediction-%x-%j.out

IDENTITY=$1

SCRIPT_DIR=$HOME/scratch/peptide-screening/predictors/PIPR
EXP_DIR=$HOME/scratch/peptide-screening/data/raw_inputs_$IDENTITY

# Load the modules
module load anaconda3 gcc

# Load the virtual env
cd $HOME
source activate piprenv
cd $SCRIPT_DIR

CUDA_VISIBLE_DEVICES=0 python3 pipr_rcnn.py $EXP_DIR/human_proteome_with_peptides.tsv $EXP_DIR/training_pipr.tsv $EXP_DIR/interactions_to_predict_pipr.tsv -load $EXP_DIR/pipr_model.model --predictions $EXP_DIR/predictions.tsv
