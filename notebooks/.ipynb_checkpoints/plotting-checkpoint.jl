### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 7c5ac541-689a-43af-8002-1f75f2658ac1
begin
	using CSV
	using DataFrames
	using Plots
	using StatsBase
	using StatsPlots
	using Pkg
	using Plots.Measures
end

# ╔═╡ a5c685cf-083c-447b-8940-4af7e8aaf069
begin
ENV["GKS_ENCODING"]="utf8"
pyplot()
end

# ╔═╡ c78ad002-b7db-11eb-3ed2-adfed9423e36
cd("/home/fcharih/Projects/PhD/peptide_screening/")

# ╔═╡ 45cdb33c-0610-4b80-87fd-b885d09a2bd1
targets_df = DataFrame(
	CSV.File("data/peptide_subset.csv")
)[!, [:compound_name, :targets]]

# ╔═╡ e7ea380f-a7b7-41a8-ae0c-9278fdf3ec19
targets = Dict()

# ╔═╡ 7b571d93-d6bc-4f91-b8da-67b462de3353
for row in eachrow(targets_df)
	targets[row["compound_name"]] = split(row["targets"], ";")
end

# ╔═╡ e136e461-d657-4f4d-9997-aa292f52ad45
targets

# ╔═╡ fdff33aa-d8ea-4a30-83ee-b0be15d6a0cd
pop!(targets, "albiglutide")

# ╔═╡ 9f280a6a-dafc-4cc9-80ea-c5d13494c25d
begin
pipr_minus = DataFrame(CSV.File("data/predictions/predictions_minus_160621_pipr.tsv", delim=" ", header=["peptide", "protein", "score"]))
pipr_plus = DataFrame(CSV.File("data/predictions/predictions_plus_160621_pipr.tsv", delim=" ", header=["peptide", "protein", "score"]))
sprint_minus = DataFrame(CSV.File("data/predictions/predictions_minus_old_sprint.tsv", delim=" ", header=["peptide", "protein","score"]))				
sprint_plus = DataFrame(CSV.File("data/predictions/predictions_plus_160621_sprint.tsv", delim=" ", header=["peptide", "protein", "score"]))
dscript = DataFrame(CSV.File("data/predictions/dscript_base_predictions.tsv", delim=" ", header=["peptide", "protein", "score"]))
end

# ╔═╡ 23782d12-f751-4bc8-a386-dac9529b44a6
scores = vcat(
	select(pipr_minus, :, :score => ByRow(x -> "pipr_minus")),
	select(pipr_plus, :, :score => ByRow(x -> "pipr_plus")),
	select(sprint_minus, :, :score => ByRow(x -> "sprint_minus")),
	select(sprint_plus, :, :score => ByRow(x -> "sprint_plus")),
	select(dscript, :, :score => ByRow(x -> "dscript")),
)

# ╔═╡ df4d9a56-38cf-4dbe-bc9c-7a1d886c95f1
 filter!(x -> x.peptide != "albiglutide", scores)

# ╔═╡ a87fc545-8e72-4dd1-a860-5914832a9b28
begin
	methods = ["dscript", "pipr_minus", "pipr_plus", "sprint_minus", "sprint_plus"]
	data = []
	for (peptide, peptide_targets) in targets
		for method in methods
			# Assign ranks
			sub = filter(x -> x.peptide == peptide && x.score_function == method, scores)
			sub["rank"] = competerank(sub[!, :score], rev = true)
			for target in peptide_targets
				selection = sub[sub.protein .== target, :]
				named_tup = Tables.rowtable(selection)
				push!(data, named_tup)
			end
		end
	end
end

# ╔═╡ f8bb23c8-fe74-47fb-8bba-670cf669b0cd
println(targets)

# ╔═╡ 750d3d9d-f86e-45e5-96c5-bd60549098b7
begin
	df2 = DataFrame(peptide=[], protein=[], score=[], score_function=[], rank=[])
	for ele in data
		element = ele[1]
		push!(df2, [element.peptide, element.protein, element.score, element.score_function, element.rank])

	end
end

# ╔═╡ e8c6f943-0226-460a-b4e0-5d0e67ff76b4
df2

# ╔═╡ bf4dc451-679e-451b-9202-e19dd7eedbbb
function plot_curve(dataframe, peptide, method, targets)
	peptide_df = filter(x -> x.peptide == peptide && x.score_function == method, dataframe)
	sort!(peptide_df, :score, rev=true)
	peptide_df["rank"] = 1:nrow(peptide_df)
	
	targets_df = filter(x -> x.protein ∈ targets[peptide], peptide_df)
	
	o2a = scatter(peptide_df.rank, peptide_df.score, framestyle=:box, legend=false, xlabel="Rank", ylabel="Score")
	scatter!(targets_df.rank, targets_df.score, color=:red, marker=8, formatter=:plain)
	
	
	o2a
end

# ╔═╡ b8d118cb-2fb8-466e-91ef-d4ea47457d6c
plot_curve(scores, "albiglutide", "sprint_minus", targets)

# ╔═╡ 57b61759-6056-493c-89f8-2eaad2347c81
ranking_matrix = vcat([
	df2[df2.score_function .== method, :]["rank"]
	for method in methods
])

# ╔═╡ 956958b7-69e3-4087-a231-b28119e8fb58
begin
	labels2 = ["PIPR (-)" "PIPR (+)" "SPRINT (-)" "SPRINT (+)" "D-SCRIPT\n(pre-trained)"]
	box_plot2 = boxplot(labels2, ranking_matrix, leg = false, framestyle=:box, outliers=false, formatter=:plain, xlabel="Method", ylabel="Rank")
	dotplot!(labels2, ranking_matrix)
end

# ╔═╡ 14dd7609-8784-4322-9200-f0fc160769b4
scores_receptors = vcat(
	select(pipr_minus, :, :score => ByRow(x -> "pipr_minus")),
	select(pipr_plus, :, :score => ByRow(x -> "pipr_plus")),
	select(sprint_minus, :, :score => ByRow(x -> "sprint_minus")),
	select(sprint_plus, :, :score => ByRow(x -> "sprint_plus")),
	select(dscript, :, :score => ByRow(x -> "dscript")),
)

# ╔═╡ fbd43006-fd8f-4eff-b6b0-ca64aa5106bd
begin
	surface_proteins = Set(readlines("data/surfaceome.txt"))
	scores_surface = filter(row -> (row.protein ∈ surface_proteins) || (row.protein ∈ targets[row.peptide]), scores)
	filter!(x -> x.peptide != "albiglutide", scores_surface)
end

# ╔═╡ c8551873-daff-42d6-9647-9ddf3b99a068


# ╔═╡ 87d71a66-9709-4f83-a68f-e34e1c0d337e
sort(filter(x -> x.peptide == "carperitide" && x.score_function == "sprint_plus", scores_surface), :score, rev=true)

# ╔═╡ b71647db-b07d-4875-9477-93ac07539cd8
length(surface_proteins)

# ╔═╡ 0d608224-268c-4386-bac4-4ab41fb47a49
filter(x -> x.peptide == "bivalirudin" && x.protein == "P00734", scores)

# ╔═╡ 3b072371-c283-45e6-a817-c5715a0583a4
begin
	methods2 = ["dscript", "pipr_minus", "pipr_plus", "sprint_minus", "sprint_plus"]
	data2 = []
	for (peptide, peptide_targets) in targets
		for method in methods2
			# Assign ranks
			sub = filter(x -> x.peptide == peptide && x.score_function == method, scores_surface)
			sub["rank"] = competerank(sub[!, :score], rev = true)
			println(peptide)
			for target in peptide_targets
				
				selection = sub[sub.protein .== target, :]
				println(selection)
				named_tup = Tables.rowtable(selection)
				push!(data2, named_tup)
			end
		end
	end
end

# ╔═╡ 37fe6fdf-ff7b-4986-84d5-820b2fc2da3e
data2

# ╔═╡ 6a0e91da-4d03-439b-b64f-5da983d2ae0c
begin
	df3 = DataFrame(peptide=[], protein=[], score=[], score_function=[], rank=[])
	for ele in data2
		
		if length(ele) == 0
			println("SKIP")
			continue
		end

		element = ele[1]

		push!(df3, [element.peptide, element.protein, element.score, element.score_function, element.rank])

	end
end

# ╔═╡ 73f0d7c2-f833-450d-a49d-6b0b63a05b54
df3

# ╔═╡ 0bafaf4b-3ea3-4cb3-88df-5331c80507e2
ranking_matrix2 = vcat([
	df3[df3.score_function .== method, :]["rank"]
	for method in methods
])

# ╔═╡ 7d9821d4-f1a3-4afc-8496-667aa67cfc9c
begin
	labels3 = ["D-SCRIPT\n(pre-trained)" "PIPR (-)" "PIPR (+)" "SPRINT (-)" "SPRINT (+)" ]
	box_plot3 = boxplot(labels3, ranking_matrix2, leg = false, framestyle=:box, outliers=false, formatter=:plain, xlabel="Method", ylabel="Rank")
	dotplot!(labels3, ranking_matrix2)
end

# ╔═╡ 45a45178-55c8-4f43-a2f6-a54309c267a0
savefig(box_plot3, "figures/surface_proteins_boxplot.pdf")

# ╔═╡ 82086054-e306-4564-87c7-5db01ee6752b
df3[df3.score_function .== "sprint_plus", :]

# ╔═╡ 592bc548-2f5b-4271-a556-a14f1f18f281
begin
	num_in_top20 = Dict()
	for method in methods
		num_in_top20[method] = nrow(filter(x -> x.rank <= 20 && x.score_function == method, df3))
	end
end

# ╔═╡ 13637bc2-ba9d-4d82-ade3-b2adae750a71
num_in_top20

# ╔═╡ ad50b95b-b876-4f97-8f41-7b84af38445c
begin
bar(["PIPR (-)", "PIPR (+)", "SPRINT (-)", "SPRINT (+)", "D-SCRIPT\n(pre-trained)"], [0, 2, 5, 15, 0], legend=false, framestyle=:box, xlabel="Method", ylabel="Count")
	savefig("figures/top20.pdf")
end

# ╔═╡ ad5e3ba8-3158-494c-bd6b-6cec1487828b


# ╔═╡ 5b7673f3-cbd1-4df9-ac22-e20c358924e0
begin
plot_curve(scores_surface, "bivalirudin", "dscript", targets)
savefig("figures/bivalirudin.pdf")
end

# ╔═╡ e7b92ae5-02e8-4e19-8594-5f929c592d19
targets

# ╔═╡ 8daa52dd-d29d-4e16-83dd-b5b7c84b4e5b


# ╔═╡ 220d4861-c3c9-409d-8487-0f60d507d674
begin
	pp = "teduglutide"
	plot_curve(scores_surface, "$pp", "sprint_plus", targets)
	savefig("figures/$(pp)_sprint_plus_surface.pdf")
end

# ╔═╡ 945781f8-c827-4cd4-8a85-cf90ae91ad12
plot_curve(scores_surface, "$pp", "dscript", targets)

# ╔═╡ 409b4844-021c-4bb2-be5a-4a6d1982cf1d
plot_curve(scores_surface, "$pp", "sprint_plus", targets)

# ╔═╡ d3c32902-d420-43f9-a179-3a91d4bd82cd
begin
	pp2 = "carperitide"
	plot_curve(scores_surface, "$(pp2)", "sprint_minus", targets)
	savefig("figures/$(pp2)_sprint_minus_surface.pdf")
	plot_curve(scores_surface, "$(pp2)", "sprint_plus", targets)
	savefig("figures/$(pp2)_sprint_plus_surface.pdf")
	plot_curve(scores_surface, "$(pp2)", "pipr_minus", targets)
	savefig("figures/$(pp2)_pipr_minus_surface.pdf")
	plot_curve(scores_surface, "$(pp2)", "pipr_plus", targets)
	savefig("figures/$(pp2)_pipr_plus_surface.pdf")
	plot_curve(scores_surface, "$(pp2)", "dscript", targets)
	savefig("figures/$(pp2)_dscript_surface.pdf")
end

# ╔═╡ b6c4aace-610b-4917-9918-7ff2730aa303
filter(x -> x.peptide == "lixisenatide", scores_surface)

# ╔═╡ 3c6371af-2f54-4a31-9746-a1f70facbd50
function get_target_ranks_for_peptide(dataframe, peptide, targets)
	df = filter(x -> x.peptide == peptide, dataframe)
	sort!(df, [:score], rev = true)
	df[!, :rank] = competerank(df[!, :score], rev = true)
	filter(row -> row.protein ∈ targets[peptide], df)[!, :rank]
end

# ╔═╡ 870c1de8-596d-454b-a661-837fb17c16a2
function get_ranks_for_method(dataframe, targets)
	ranks = []
	for peptide in keys(targets)
		peptide_ranks = get_target_ranks_for_peptide(dataframe, peptide, targets)
		append!(ranks, peptide_ranks)
	end
	return ranks
end

# ╔═╡ 54d7153a-3a86-4f5f-92f7-ac9b48c9eca1
get_ranks_for_method(dscript, targets)

# ╔═╡ f6ba3f69-fee5-42b7-bbc9-3113da34bd9a
ranks = DataFrame([
	get_ranks_for_method(pipr_minus, targets),
	get_ranks_for_method(pipr_plus, targets),
	get_ranks_for_method(sprint_minus, targets),
	get_ranks_for_method(sprint_plus, targets),
	get_ranks_for_method(dscript, targets)
])

# ╔═╡ 35dc5c45-6d0e-4b73-915d-56f2c8e06795
rename!(ranks, Symbol.(["pipr_minus", "pipr_plus", "sprint_minus", "sprint_plus", "dscript"]))

# ╔═╡ 8f7bddcd-71a5-4c3b-ac27-36c40c429455
begin
	ranks_matrix = Matrix{Int64}(ranks)
	labels = ["PIPR (-)" "PIPR (+)" "SPRINT (-)" "SPRINT (+)" "D-SCRIPT\n(pre-trained)"]
	box_plot = boxplot(labels, ranks_matrix, leg = false, framestyle=:box, outliers=false, formatter=:plain, xlabel="Method", ylabel="Rank")
	dotplot!(labels, ranks_matrix)
end

# ╔═╡ bd2e27af-efdd-4512-a262-c705d925666c
savefig(box_plot, "figures/boxplot.pdf")

# ╔═╡ 0124ff9f-cc00-4dd4-bd34-c9fe7df1c2fd
function get_hist_series(ranks)
	x = range(1, 20400, step=1)
	y = []
	for i in x
		push!(y, length(filter(x -> x < i, ranks)))
	end
	return (x, y)
end

# ╔═╡ 2f6e1830-1b02-458e-9ae1-38908222917a
begin
rising_curve = plot(get_hist_series(ranks.pipr_minus), label="PIPR (-)", framestyle=:box, xlabel="Rank cutoff", ylabel="Targets below cutoff", formatter=:plain)
plot!(get_hist_series(ranks.pipr_plus), label="PIPR (+)")
plot!(get_hist_series(ranks.sprint_minus), label="SPRINT (-)")
plot!(get_hist_series(ranks.sprint_plus), label="SPRINT (+)")
plot!(get_hist_series(ranks.dscript), label="D-SCRIPT (pre-trained)")
	
plot!(
    [
	get_hist_series(ranks.pipr_minus),
	get_hist_series(ranks.pipr_plus),
	get_hist_series(ranks.sprint_minus),
	get_hist_series(ranks.sprint_plus),
	get_hist_series(ranks.dscript)
	],
    inset = bbox(0.3, 0.2, 0.3, 0.3, :center),
    subplot = 2,
	framestyle=:box,
	formatter=:plain,
	legend=nothing,
	xrange=(0,50),
	yrange=(0,20)
)
end

# ╔═╡ 6cb631f5-bcc1-45f6-b974-e7015b888384
savefig(rising_curve, "figures/targets_below_rank_plot.pdf")

# ╔═╡ 60b2faaa-71d2-4f77-850a-4910c72e04d3
histogram(sprint_plus.score, bins=100, yrange=(0,1000))

# ╔═╡ Cell order:
# ╠═7c5ac541-689a-43af-8002-1f75f2658ac1
# ╠═a5c685cf-083c-447b-8940-4af7e8aaf069
# ╠═c78ad002-b7db-11eb-3ed2-adfed9423e36
# ╠═45cdb33c-0610-4b80-87fd-b885d09a2bd1
# ╠═e7ea380f-a7b7-41a8-ae0c-9278fdf3ec19
# ╠═7b571d93-d6bc-4f91-b8da-67b462de3353
# ╠═e136e461-d657-4f4d-9997-aa292f52ad45
# ╠═fdff33aa-d8ea-4a30-83ee-b0be15d6a0cd
# ╠═9f280a6a-dafc-4cc9-80ea-c5d13494c25d
# ╠═23782d12-f751-4bc8-a386-dac9529b44a6
# ╠═df4d9a56-38cf-4dbe-bc9c-7a1d886c95f1
# ╠═a87fc545-8e72-4dd1-a860-5914832a9b28
# ╠═f8bb23c8-fe74-47fb-8bba-670cf669b0cd
# ╠═750d3d9d-f86e-45e5-96c5-bd60549098b7
# ╠═e8c6f943-0226-460a-b4e0-5d0e67ff76b4
# ╠═bf4dc451-679e-451b-9202-e19dd7eedbbb
# ╠═b8d118cb-2fb8-466e-91ef-d4ea47457d6c
# ╠═57b61759-6056-493c-89f8-2eaad2347c81
# ╠═956958b7-69e3-4087-a231-b28119e8fb58
# ╠═14dd7609-8784-4322-9200-f0fc160769b4
# ╠═fbd43006-fd8f-4eff-b6b0-ca64aa5106bd
# ╠═c8551873-daff-42d6-9647-9ddf3b99a068
# ╠═87d71a66-9709-4f83-a68f-e34e1c0d337e
# ╠═945781f8-c827-4cd4-8a85-cf90ae91ad12
# ╠═b71647db-b07d-4875-9477-93ac07539cd8
# ╠═0d608224-268c-4386-bac4-4ab41fb47a49
# ╠═3b072371-c283-45e6-a817-c5715a0583a4
# ╠═37fe6fdf-ff7b-4986-84d5-820b2fc2da3e
# ╠═6a0e91da-4d03-439b-b64f-5da983d2ae0c
# ╠═73f0d7c2-f833-450d-a49d-6b0b63a05b54
# ╠═0bafaf4b-3ea3-4cb3-88df-5331c80507e2
# ╠═7d9821d4-f1a3-4afc-8496-667aa67cfc9c
# ╠═45a45178-55c8-4f43-a2f6-a54309c267a0
# ╠═82086054-e306-4564-87c7-5db01ee6752b
# ╠═592bc548-2f5b-4271-a556-a14f1f18f281
# ╠═13637bc2-ba9d-4d82-ade3-b2adae750a71
# ╠═ad50b95b-b876-4f97-8f41-7b84af38445c
# ╠═ad5e3ba8-3158-494c-bd6b-6cec1487828b
# ╠═5b7673f3-cbd1-4df9-ac22-e20c358924e0
# ╠═e7b92ae5-02e8-4e19-8594-5f929c592d19
# ╠═8daa52dd-d29d-4e16-83dd-b5b7c84b4e5b
# ╠═220d4861-c3c9-409d-8487-0f60d507d674
# ╠═409b4844-021c-4bb2-be5a-4a6d1982cf1d
# ╠═d3c32902-d420-43f9-a179-3a91d4bd82cd
# ╠═b6c4aace-610b-4917-9918-7ff2730aa303
# ╠═3c6371af-2f54-4a31-9746-a1f70facbd50
# ╠═870c1de8-596d-454b-a661-837fb17c16a2
# ╠═54d7153a-3a86-4f5f-92f7-ac9b48c9eca1
# ╠═f6ba3f69-fee5-42b7-bbc9-3113da34bd9a
# ╠═35dc5c45-6d0e-4b73-915d-56f2c8e06795
# ╠═8f7bddcd-71a5-4c3b-ac27-36c40c429455
# ╠═bd2e27af-efdd-4512-a262-c705d925666c
# ╠═0124ff9f-cc00-4dd4-bd34-c9fe7df1c2fd
# ╠═2f6e1830-1b02-458e-9ae1-38908222917a
# ╠═6cb631f5-bcc1-45f6-b974-e7015b888384
# ╠═60b2faaa-71d2-4f77-850a-4910c72e04d3
