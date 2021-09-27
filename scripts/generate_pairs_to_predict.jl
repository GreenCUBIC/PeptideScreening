using FASTX
using Random
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--proteins", "-p"
        required = true
        arg_type = String
        help = "Path to the FASTA/tsv file containing the protein sequences."
    "--peptides", "-e"
        required = true
        arg_type = String
        help = "Path to the FASTA/tsv file containing the peptide sequences."
    "--output", "-o"
        required = true
        arg_type = String
        help = "Path to the tsv file listing the interactions to score for PIPR."
end

args = parse_args(ARGS, s)

PROTEINS = args["proteins"]
PEPTIDES = args["peptides"]
OUTPUT = args["output"]
PROTEIN_FASTA = occursin("fasta", PROTEINS) ? true : false
PEPTIDES_FASTA = occursin("fasta", PEPTIDES) ? true : false

if PROTEIN_FASTA
    protein_ids = [FASTA.identifier(record) for record in FASTA.Reader(open(PROTEINS))]
else
    protein_ids = [split(line, "\t")[1] for line in readlines(PROTEINS)]
end

if PEPTIDES_FASTA
    peptide_ids = [FASTA.identifier(record) for record in FASTA.Reader(open(PEPTIDES))]
else
    peptide_ids = [split(line, "\t")[1] for line in readlines(PEPTIDES)]
end

open(OUTPUT, "w") do io
    write(io, "p1\tp2\tlabel\n")
	for peptide in peptide_ids
            for protein in protein_ids
		write(io, "$peptide\t$protein\t$(rand([0,1]))\n")
            end
	end
end
