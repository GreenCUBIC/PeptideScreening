using FASTX
using DataFrames
using CSV
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--output", "-o"
        help = "Path to the tsv file containing a list of all the human interactions (conservative) from the BIOGRID database."
        arg_type = String
        required = true
    "--biogrid_version", "-b"
        help = "Version of BIOGRID to download (default: 4.3.196)"
        arg_type = String
        required = false
        default = "4.3.196"
end

args = parse_args(ARGS, s)

DESTINATION = args["output"]
BIOGRID_VERSION = args["biogrid_version"]
BIOGRID_URL = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-$BIOGRID_VERSION/BIOGRID-ALL-$BIOGRID_VERSION.tab3.zip"

# Download the interactions only if not already available
run(`curl $BIOGRID_URL -o biogrid.tab3.zip`)

# Unzip
run(`unzip biogrid.tab3.zip`)
run(`rm biogrid.tab3.zip`)

#setenv(`unzip biogrid.tab3.zip`, dir="/tmp")

# Load the interaction data
interactions_df = DataFrame(CSV.File("BIOGRID-ALL-$BIOGRID_VERSION.tab3.txt"))
run(`rm BIOGRID-ALL-$BIOGRID_VERSION.tab3.txt`)

# Filter out non-human
filter!(row -> row["Organism ID Interactor A"] == 9606 && row["Organism ID Interactor B"] == 9606, interactions_df)

# Only keep the interactions discovered with stringent methods (high-quality interactions)
CONSERVATIVE_SYSTEMS = Set([
      "Two-hybrid",
      "Affinity Capture-MS",
      "Affinity Capture-Western",
      "Reconstituted Complex",
      "Affinity Capture-Luminescence",
      "Co-crystal Structure",
      "Far Western",
      "FRET",
      "Protein-peptide",
      "Co-localization",
      "Affinity Capture-RNA",
      "Co-purification"
])

# Only keep interactions detected with one of the conservative systems
filter!(row -> row["Experimental System"] ∈ CONSERVATIVE_SYSTEMS, interactions_df)

# Get rid of useless columns
interactions = select(
    interactions_df[:, [
        "SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B", "Publication Source"
    ]],
    "SWISS-PROT Accessions Interactor A" => "prot_a",
    "SWISS-PROT Accessions Interactor B" => "prot_b",
    "Publication Source" => "source"
)

# Count the number of times the interactions have been validated in different sources
sources = Dict()
pair_counts = Dict()

for pair in eachrow(interactions)
    p1 = min(pair["prot_a"], pair["prot_b"])
    p2 = max(pair["prot_a"], pair["prot_b"])
    source = pair["source"]
    
    pair_tuple = (p1, p2)
    
    if pair_tuple ∉ keys(pair_counts)
        pair_counts[pair_tuple] = 1
        sources[pair_tuple] = [source]
    elseif pair_tuple ∈ keys(pair_counts) && source ∉ sources[pair_tuple]
        pair_counts[pair_tuple] += 1
        push!(sources[pair_tuple], source)
    else
        continue
    end
end

# Only keep pairs that are present more than once (multiple lines of evidence)
duplicated_pairs = [k for (k, v) in pair_counts if v >= 2]

open(DESTINATION, "w") do io
    write(io, join(["$(pair[1])\t$(pair[2])" for pair in duplicated_pairs if pair[1] != "-" && pair[2] != "-"], "\n"))
end
