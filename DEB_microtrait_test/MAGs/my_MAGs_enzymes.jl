using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
df_MAGs         = CSV.read("files/MAGs2traits_DEBinput.csv", DataFrame, missingstring="N/A")
n_MAGs          = nrow(df_MAGs)
println("Loaded $n_MAGs MAGs")
########################################

########################################
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_MAGs.Genome_size))
zh              = df_MAGs.z_hydrolases./df_MAGs.Genome_size*1e6
alpha_X         = 1e-2*(df_MAGs.z_hydrolases./df_MAGs.Genome_size*1e6)./maximum(df_MAGs.z_hydrolases./df_MAGs.Genome_size*1e6)
println("Median hydrolase density: $(median(zh))")
########################################

########################################
save("files/my_MAGs_enzymes.jld", "zh", zh, "alpha", alpha_X)
println("Saved: files/my_MAGs_enzymes.jld")
########################################
