using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

########################################
# I/O
df_isolates     = CSV.read("files/isolates2traits_DEBinput.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("files/exudation_properties.csv", DataFrame, missingstring="N/A")
df_metabolites.Formula = convert.(String, df_metabolites.Formula)
n_isolates      = nrow(df_isolates)
println("Loaded $n_isolates isolates")

# Glucose only — index 29 in exudation_properties.csv
glucose_idx     = 29
println("Running assimilation for glucose only (index $glucose_idx: $(df_metabolites.Name[glucose_idx]))")
########################################

########################################
# isolate traits
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)

z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6), 1, n_isolates)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6), 1, n_isolates)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1, n_isolates)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6), 1, n_isolates)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6), 1, n_isolates)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6), 1, n_isolates)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

y_EM            = ones(size(V_cell,1))
########################################

########################################
# calc transporter density for glucose only
rho_ps          = zeros(1, size(V_cell,1))
y_DEs           = zeros(1, size(V_cell,1))
failed_isolates = Int[]

j = glucose_idx
for i in 1:size(V_cell,1)
    try
        find_rho(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
        rho_p = Roots.find_zero(find_rho, 1.0)
        closure = genome_distr[:,i]./sum(genome_distr[:,i])
        rho_ps[1,i] = rho_p[1].*closure[1]
        y_DE = DEBmicroTrait.yield_transporter_density_cost(rho_ps[1,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
        y_DEs[1,i] = y_DE[1]
    catch e
        rho_ps[1,i] = 1e-12
        y_DEs[1,i] = 0.0
        push!(failed_isolates, i)
    end

    if i % 500 == 0
        println("Processing isolate $i of $n_isolates... ($(length(failed_isolates)) failed so far)")
        save("files/my_isolates_assimilation_glucose_partial.jld", "rho", rho_ps, "yDE", y_DEs, "progress", i)
    end
end

println("\nCompleted all $n_isolates isolates")
println("Total convergence failures: $(length(failed_isolates))")
if length(failed_isolates) > 0
    println("Failed isolate indices: $failed_isolates")
end

# Clean up any remaining problematic values
rho_ps[rho_ps.==0.0] .= 1e-12
rho_ps[isnan.(rho_ps)] .= 1e-12
rho_ps[isinf.(rho_ps)] .= 1e-12
rho_ps[rho_ps.<0.0] .= 1e-12
y_DEs[isnan.(y_DEs)] .= 0.0
y_DEs[isinf.(y_DEs)] .= 0.0
println("Median transporter density: $(median(rho_ps))")

# N_C for glucose only
N_C = zeros(1)
N_C[1] = DEBmicroTrait.extract_composition(df_metabolites.Formula[glucose_idx])[1]
println("N_C for glucose: $(N_C[1])")

N_SB  = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, rho_ps, Min_gen_time, Gram_stain)
Vmax  = @. 180.0*60^2*N_SB.*N_C
D_S   = DEBmicroTrait.aqueous_diffusivity([df_metabolites.Molecular_weight[glucose_idx]])
K_D   = DEBmicroTrait.specific_reference_affinity(V_cell, rho_ps, D_S)
a_s   = Vmax./K_D

########################################
# I/O — final save
save("files/my_isolates_assimilation_glucose.jld", "rho", rho_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C, "failed", failed_isolates)
println("Saved: files/my_isolates_assimilation_glucose.jld")
########################################
