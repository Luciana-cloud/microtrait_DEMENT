using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

########################################
# I/O
df_isolates     = CSV.read("files/isolates2traits_DEBinput.csv", DataFrame, missingstring="N/A")
n_isolates      = nrow(df_isolates)
println("Loaded $n_isolates isolates")
########################################

########################################
Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
########################################

########################################
k_M             = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell, Min_gen_time, Gram_stain)
k_M_med         = median(k_M)
println("Median maintenance rate: $k_M_med")
########################################

########################################
# I/O
save("files/my_isolates_maintenance.jld", "kM", k_M)
println("Saved: files/my_isolates_maintenance.jld")
########################################
