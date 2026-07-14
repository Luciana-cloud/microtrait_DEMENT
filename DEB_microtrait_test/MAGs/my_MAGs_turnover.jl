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
Genome_size     = convert(Array{Float64,1}, df_MAGs.Genome_size)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time    = df_MAGs.Min_gen_time
gmax            = log(2)./Min_gen_time
Gram_stain      = convert(Array{String,1}, df_MAGs.gram_stain)
########################################

########################################
gamma_V_0       = DEBmicroTrait.max_specific_death_rate(Min_gen_time::Vector{Float64})
println("Median max death rate: $(median(gamma_V_0))")
########################################

########################################
dry_mass        = 0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
rho_bulk        = 1.21
Bio_0           = 1e9*1e6*rho_bulk*dry_mass./12.011
gamma_V_1       = median(Bio_0)
########################################

########################################
save("files/my_MAGs_turnover.jld", "gV0", gamma_V_0, "gV1", gamma_V_1)
println("Saved: files/my_MAGs_turnover.jld")
########################################
