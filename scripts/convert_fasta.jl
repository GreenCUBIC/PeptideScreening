#! /usr/local/bin/julia
using ArgParse
using FASTX

s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        arg_type = String
        required = true
        help = "Path to the input file (FASTA or tsv)."
    "--output", "-o"
        arg_type = String
        required = true
        help = "Path to the output file."
end

args = parse_args(ARGS, s)

IS_FASTA = occursin("fasta", args["input"])

if IS_FASTA
    records = [rec for rec in FASTA.Reader(open(args["input"]))]
    fasta_content = ["$(FASTA.identifier(rec))\t$(FASTA.sequence(rec))" for rec in records]
    open(args["output"], "w") do io
        write(io, join(fasta_content, "\n"))
    end
else
    proteins = [split(line, "\t") for line in readlines(args["input"])]
    open(args["output"], "w") do io
        for protein in proteins
            write(io, ">$(protein[1])\n$(protein[2])\n")
        end
    end
end
