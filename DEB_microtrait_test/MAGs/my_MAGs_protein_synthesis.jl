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
rrn_copies      = convert(Array{Float64,1}, df_MAGs.rRNA_genes)
Min_gen_time    = df_MAGs.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_MAGs.gram_stain)
########################################

########################################
gmax            = log(2)./Min_gen_time
V_p             = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_r             = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
k_E             = DEBmicroTrait.translation_power(V_p, V_r, Min_gen_time)
println("Median translation power: $(median(k_E))")
########################################

########################################
y_EV            = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
println("Median translation efficiency: $(median(y_EV))")
########################################

########################################
save("files/my_MAGs_protein_synthesis.jld", "kE", k_E, "yEV", y_EV, "mingt", Min_gen_time)
println("Saved: files/my_MAGs_protein_synthesis.jld")
########################################
