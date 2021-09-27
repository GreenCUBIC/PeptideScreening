using Plots
using Plots.Measures
using CSV
using DataFrames
using ArgParse

pyplot()

s = ArgParseSettings()
@add_arg_table s begin
    "--scores", "-s"
        required = true
        arg_type = String
        help = "Path to the file containing the scores in tsv format (peptide protein score)."
    "--peptides", "-p"
        required = true
        arg_type = String
        help = "Path to the csv file containing the peptides data."
    "--output", "-o"
        required = true
        arg_type = String
        help = "Output file where the plots are to be saved."
end

args = parse_args(ARGS, s)

PEPTIDES_FILEPATH = args["peptides"]

SCORES_FILEPATH = args["scores"]

function load_peptide_targets(peptide_filepath::String)
	df = DataFrame(CSV.File(peptide_filepath))

	peptides = Dict{String, Array{String}}()
	for row in eachrow(df)
		peptides[row["compound_name"]] = split(row["targets"], ";")
	end

	peptides
end

function generate_curves(scores_filepath::String)
	df = DataFrame(CSV.File(scores_filepath, header=false, delim=" "))
	peptide_groups = groupby(df, :Column1)
	peptide_curves = Dict()
    println(names(peptide_groups))
	for group in peptide_groups
		peptide = group[1,:]["Column1"]
		sorted = sort(group, :Column3, rev=true)
		sorted[!, :rank] = 1:nrow(sorted)
		peptide_curves[peptide] = sorted[!,[:rank,:Column2,:Column3]]
	end
	peptide_curves
end

function plot_curve(peptide::String, curves::Dict{Any, Any}, targets::Dict{String, Array{String}})
	df = curves[peptide]
	non_targets = filter(x -> x["Column2"] ∉ targets[peptide], df)
	targets = filter(x -> x["Column2"] ∈ targets[peptide], df)
	scatter(
		non_targets[!, :rank],
		non_targets[!, :Column3],
	  	title="O2A Curve ($peptide)",
		label="Non-target(s)",
		xlabel="Rank", ylabel="Score",
        margin=10mm
       )
	scatter!(
		targets[!, :rank],
		targets[!, :Column3],
		markersize=6,
		color="red", label="Target(s)")
end

println("Loading the peptide data...")
targets = load_peptide_targets(PEPTIDES_FILEPATH)

println("Generating the O2A curves...")
curves = generate_curves(SCORES_FILEPATH)

println("Plotting the curves...")
plots = [plot_curve(peptide, curves, targets) for peptide in keys(targets)]
plot(plots..., layout = length(plots), legend = false, size = (2000,2000), margin=20mm)

println("Saving the plots to $(args["output"])...")
savefig(args["output"])
