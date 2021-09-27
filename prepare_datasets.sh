#! /bin/bash

# Create data directory
mkdir -p data

# Download the proteome (in fasta)
julia scripts/download_human_sequences.jl -o data/human_proteome.fasta

# Append the peptides
cat data/human_proteome.fasta data/peptide_sequences.fasta > data/proteome_and_peptide_sequences.fasta

# Make a TSV copy
julia scripts/convert_fasta.jl -i data/proteome_and_peptide_sequences.fasta -o data/proteome_and_peptide_sequences.tsv

# Download the conservative interactions (as per Dick et al. 2018)
julia scripts/download_conservative_interactions.jl -o data/conservative_interactions.tsv

# Identify the endogenous equivalents of the peptides
julia scripts/identify_endogenous_equivalents.jl -s data/human_proteome.fasta -p data/peptide_sequences.tsv -o data/endogenous_equivalents.txt

# Create the different training data
julia scripts/split_pairs.jl -s data/human_proteome.fasta -i data/conservative_interactions.tsv -o data -d data/endogenous_equivalents.txt

# Generate a list of pairs for which the scores are sought (for prediction)
julia scripts/generate_pairs_to_predict.jl -p data/human_proteome.fasta -e data/peptide_sequences.tsv -o data/interactions_to_predict_pipr.tsv
awk 'NF{NF-=1};1' data/interactions_to_predict_pipr.tsv | tail -n +2 > data/interactions_to_predict_dscript.tsv # Remove header and random class
