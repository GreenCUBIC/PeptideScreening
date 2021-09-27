using FASTX

DATA_DIR = "../data"
SEQ_DIR = "$DATA_DIR/sequences"
INT_DIR = "$DATA_DIR/interactions"
PIPR_DIR = "$DATA_DIR/pipr"

TRAINING_SET_PATH = "$INT_DIR/training.tsv"
TEST_SET_PATH = "$INT_DIR/test.tsv"
SEQUENCES_PATH = "$SEQ_DIR/human_sequences_clustered.fasta"
ALL_SEQUENCES_PATH = "$SEQ_DIR/human_sequences.fasta"
PEPTIDES_PATH = "$SEQ_DIR/peptide_sequences.tsv"
FORMATTED_SEQUENCE_PATH = "$PIPR_DIR/formatted_sequences_with_peptides.tsv"

println("Load the protein sequences...")
proteins = []
for record in FASTA.Reader(open(SEQUENCES_PATH, "r"))
    protein_id = FASTA.identifier(record)
    push!(proteins, "$(protein_id)\t$(FASTA.sequence(record))")
end

println("Loading the peptide sequences...")
peptides = [split(p, "\t") for p in readlines(PEPTIDES_PATH)]
all_proteins = [proteins; ["$(peptide[1])\t$(peptide[2])" for peptide in peptides]]

println("Saving the sequences...")
formatted_sequences = join(all_proteins, "\n")
write(open(FORMATTED_SEQUENCE_PATH, "w"), formatted_sequences)

human_sequences_path = "$PIPR_DIR/human_sequences.tsv"
cp(ALL_SEQUENCES_PATH, human_sequences_path, force=true)

train_target_path = "$PIPR_DIR/train.tsv"
cp(TRAINING_SET_PATH, train_target_path, force=true)

test_target_path = "$PIPR_DIR/test.tsv"
cp(TEST_SET_PATH, test_target_path, force=true)

peptide_interactions = "$PIPR_DIR/peptide_interactions.tsv"
open(peptide_interactions, "w") do io
    write(io, "p1\tp2\tlabel\n")
    for peptide in peptides
        for protein in proteins
            peptide_id = peptide[1]
            protein_id = split(protein, "\t")[1]
            random_class = rand([0,1])
            write(io, "$(peptide_id)\t$(protein_id)\t$(random_class)\n")
        end
    end
end
