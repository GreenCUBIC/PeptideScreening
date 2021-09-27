#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --gpus-per-node=4
#SBATCH --nodes=1
#SBATCH -p compute_full_node
#SBATCH --account=soscip-3-091
#SBATCH --output=/home/j/jrgreen/fcharih/scratch/peptide-screening/logs/train_dscript-%x-%j.out

IDENTITY=$1

EXP_DIR=$HOME/scratch/peptide-screening/data/raw_inputs_$IDENTITY
EMBEDDINGS=$HOME/scratch/peptide-screening/data/outputs/human_proteome_with_peptides.embed

# Load the modules
module load anaconda3 gcc cuda cudnn

# Load the virtual env
cd $HOME
source activate dscriptenv

dscript train --train $EXP_DIR/training_dscript.tsv --val $EXP_DIR/validation_dscript.tsv --embedding $EMBEDDINGS --save-prefix $EXP_DIR/dscript_model.model --device 0 --batch-size 4
