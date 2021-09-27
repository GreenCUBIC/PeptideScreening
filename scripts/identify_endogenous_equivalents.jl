using FASTX
using BioTools.BLAST
using BioSequences
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--sequences", "-s"
        arg_type = String
        help = "Path to the FASTA file with the human sequences."
        required = true
    "--peptides", "-p"
        help = "Path to the tsv file with the peptides sequences."
        arg_type = String
        required = true
    "--output", "-o"
        help = "Path to the text file with generate with one endogenous peptide equivalent per line."
        arg_type = String
        required = true
end

args = parse_args(ARGS, s)

# Retrive the peptides
peptides = [split(line, "\t") for line in readlines(args["peptides"])]

println("Blasting them peptides against the proteome...")
proteins_to_exclude = []
for peptide in peptides
    println("BLASTING $(peptide[1])...")
    seq = BioSequences.AminoAcidSequence(peptide[2])
    results = blastp(seq, args["sequences"])
    println("Found $(length(results)) results...")
    for result in results
        println(result.hitname)
        push!(proteins_to_exclude, result.hitname)
    end
end

println("Saving the list of endogenous peptides equivalent to the FDA-approved peptides to a file.")
open(args["output"], "w") do io
    write(io, join(Set(proteins_to_exclude), "\n"))
end


