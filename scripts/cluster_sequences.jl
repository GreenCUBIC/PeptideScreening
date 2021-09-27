using FASTX
using Random
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--sequences", "-s"
        arg_type = String
        required = true
        help = "Path to the FASTA file containing the protein sequences."
    "--interactions", "-i"
        arg_type = String
        required = true
        help = "Path to the tsv file containing the PPIs."
    "--output", "-o"
        arg_type = String
        help = "Path to the FASTA output file containing the reduced set of sequences."
        required = true
        arg_type = String
    "--threshold", "-t"
        default = 0.70
        arg_type = Float64
        help = "The identity threshold to be used for clustering"
        required = false
end

args = parse_args(ARGS, s)

function compute_word_size(threshold::Float64)
    if threshold >= 0.7
        return 5
    elseif threshold > 0.6
        return 4
    elseif threshold > 0.5
        return 3
    else
        return 2
    end
end

SIMILARITY_THRESHOLD = (args["threshold"], compute_word_size(args["threshold"])) # threshold, word size

# Cluster the sequences
println("Clustering the sequences...")
cmd = `cd-hit -i $(args["sequences"]) -o /tmp/temp_clusters.fasta -n $(SIMILARITY_THRESHOLD[2]) -c $(SIMILARITY_THRESHOLD[1]) -T 8`
run(cmd)

# Interaction counting
println("Counting the number of interactions proteins are involved in...")

interaction_counts = Dict()
for pair in readlines(args["interactions"])
    split_pair = split(pair, "\t")
    for protein in split_pair
        if !haskey(interaction_counts, protein)
            interaction_counts[protein] = 1
        else
            interaction_counts[protein] += 1
        end
    end
end

# Add a value of 0 to proteins which are involed in 0 interactions
for protein in FASTA.Reader(open(args["sequences"], "r"))
	id = FASTA.identifier(protein)
	if !haskey(interaction_counts, id)
		interaction_counts[id] = 0
	end
end

println("Identifying the cluster representatives based on interaction counts...")
representatives = []
cluster_file = read("/tmp/temp_clusters.fasta.clstr", String)
clusters = split(cluster_file, ">C")[2:end]

# For every cluster, identify the proteins in the cluster
for cluster in clusters
    lines = split(cluster, "\n")
    if length(lines) == 4
        representative = lines[2]
        push!(representatives, representative)
    else
        best_candidate = ""
        best_score = -1
        for candidate in lines
            m = match(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", candidate)
            if m == nothing
                continue
            else
                if interaction_counts[m.match] > best_score
                    best_candidate = m.match
                    best_score = interaction_counts[m.match]
                end
            end
        end
        push!(representatives, best_candidate)
    end
end

println("Saving the cluster representatives...")
representative_set = Set(representatives)
selected_proteins = []
for protein in FASTA.Reader(open(args["sequences"], "r"))
    id = FASTA.identifier(protein)
    if in(id, representative_set)
        push!(selected_proteins, protein)
    end
end

w = open(FASTA.Writer, args["output"])
for p in selected_proteins
    write(w, p)
end
