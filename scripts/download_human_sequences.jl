using FASTX
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--destination", "-o"
        arg_type = String
        required = true
        help = "Destination of the FASTA file."
end
args = parse_args(ARGS, s)

# Globals
DESTINATION = args["destination"]
DIRECTORY = abspath(dirname(DESTINATION))
SWISSPROT_ENDPOINT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

# Dowload the dataset
file_destination = "$DIRECTORY/uniprot_sprot.fasta.gz"
run(`curl $SWISSPROT_ENDPOINT -o $file_destination`)

# Unzip
println("Unziping the sequences...")
run(setenv(`gzip -f -d uniprot_sprot.fasta.gz`, dir=DIRECTORY))

# Filter the sequences
println("Filtering out non-human sequences...")
human_sequences = [record for record in FASTA.Reader(open("$DIRECTORY/uniprot_sprot.fasta")) if occursin("OX=9606", FASTA.description(record))]

# Dump to a file
println("Dumping sequences to $DESTINATION...")
open(DESTINATION, "w") do io
    fasta_content = join([">$(split(FASTA.identifier(seq), "|")[2])\n$(FASTA.sequence(seq))" for seq in human_sequences], "\n")
    write(io, fasta_content * "\n")
end
