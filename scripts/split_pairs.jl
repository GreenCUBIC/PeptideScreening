using FASTX
using Random
using DataFrames
using ProgressMeter
using CSV
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--sequences", "-s"
        required = true
        arg_type = String
        help = "Path to the FASTA file containing the protein sequences."
    "--interactions", "-i"
        required = true
        arg_type = String
        help = "Path to the tsv file containing the PPIs."
    "--output_dir", "-o"
        required = true
        help = "Path to the directory where the training and test pairs are to be saved."
        arg_type = String
    "--random_seed", "-e"
        required = false
        help = "The random seed to use to shuffle the interactions."
        arg_type = Int
        default = 42
    "--neg_to_pos", "-r"
        required = false
        help = "Ratio of negatives to positives to put in the training sets (r value)."
        arg_type = Float64
        default = 1.0
    "--remove", "-d"
        required = true
        help = "Path to a list of proteins (.txt) which cannot be involved in an interaction."
        arg_type = String
end

function generate_random_pairs(k::Int64, protein_ids::Set{String}, exclude::Any)
    random_pairs = Set()

    p = Progress(k, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:green)
    while length(random_pairs) < k
        p1 = rand(protein_ids)
        p2 = rand(protein_ids)
        if !in([p1, p2], exclude) && !in([p2, p1], exclude) && !in([p1, p2], random_pairs) && !in([p2, p1], random_pairs)
            push!(random_pairs, [p1, p2, 0])
        end
        next!(p)
    end

    return random_pairs
end

function write_pairs_to_file_sprint(filepath::String, pairs::Any)
    open(filepath, "w") do io
        for pair in pairs
            write(io, "$(pair[1])\t$(pair[2])\n")
        end
    end
end


function write_pairs_to_file(filepath::String, pairs::Any, header::Bool)
    open(filepath, "w") do io
    if header
        write(io, "p1\tp2\tlabel\n")
    end
    for pair in pairs
        write(io, "$(pair[1])\t$(pair[2])\t$(pair[3])\n")
    end
    end
end

function format_pos_pair(pair)
    [pair[1], pair[2], 1]
end

args = parse_args(ARGS, s)

PROTEIN_SEQUENCES = args["sequences"]
INTERACTIONS = args["interactions"]
REMOVE = args["remove"]
RANDOM_SEED = args["random_seed"]
OUTPUT_DIR = args["output_dir"]
NEG_TO_POS = args["neg_to_pos"]

# Load the protein sequences and the interactions for which the
# sequences are available
proteins = [protein for protein in FASTA.Reader(open(PROTEIN_SEQUENCES, "r"))]
protein_ids = Set([FASTA.identifier(p) for p in proteins])
pairs = [format_pos_pair(split(line, "\t")) for line in readlines(INTERACTIONS)]

# Must NOT involve an endogenous peptide equivalent
endogenous_equivalents = Set(readlines(REMOVE))
pairs_without_endogenous = [[pair[1], pair[2], 1] for pair in pairs if pair[1] ∉ endogenous_equivalents && pair[2] ∉ endogenous_equivalents]
println("There are $(length(pairs_without_endogenous))/$(length(pairs)) not involving an equivalent...")

# Create training data for SPRINT
write_pairs_to_file_sprint("$OUTPUT_DIR/sprint_training_plus.tsv", pairs)
write_pairs_to_file_sprint("$OUTPUT_DIR/sprint_training_minus.tsv", pairs_without_endogenous)

# Create data for PIPR
pairs_to_exclude = Set([split(line, "\t") for line in readlines(INTERACTIONS)])
pipr_negatives_plus = generate_random_pairs(length(pairs), protein_ids, pairs_to_exclude)
pipr_negatives_minus = generate_random_pairs(length(pairs_without_endogenous), protein_ids, pairs_to_exclude)
pipr_plus = union(Set(pairs), pipr_negatives_plus)
pipr_minus = union(Set(pairs_without_endogenous), pipr_negatives_minus)

write_pairs_to_file("$OUTPUT_DIR/pipr_training_plus.tsv", pipr_plus, true)
write_pairs_to_file("$OUTPUT_DIR/pipr_training_minus.tsv", pipr_minus, true)
write_pairs_to_file("$OUTPUT_DIR/pipr_test_bogus.tsv", [["corticotropin", "corticotropin", 1]], true)
