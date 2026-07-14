using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

########################################
# I/O
df_MAGs         = CSV.read("files/MAGs2traits_DEBinput.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("files/exudation_properties.csv", DataFrame, missingstring="N/A")
df_metabolites.Formula = convert.(String, df_metabolites.Formula)
n_MAGs          = nrow(df_MAGs)
println("Loaded $n_MAGs MAGs")

# Glucose only — index 29 in exudation_properties.csv
glucose_idx     = 29
println("Running assimilation for glucose only (index $glucose_idx: $(df_metabolites.Name[glucose_idx]))")
# println("Starting at: $(Dates.now())")
########################################

########################################
# MAG traits
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_MAGs.Genome_size))
Min_gen_time    = df_MAGs.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_MAGs.gram_stain)
rrn_copies      = convert(Array{Float64,1}, df_MAGs.rRNA_genes)

z_sugars        = reshape(convert(Array{Float64,1}, df_MAGs.z_sugars./df_MAGs.Genome_size*1e6), 1, n_MAGs)
z_organics      = reshape(convert(Array{Float64,1}, df_MAGs.z_organic_acids./df_MAGs.Genome_size*1e6), 1, n_MAGs)
z_aminos        = reshape(convert(Array{Float64,1}, df_MAGs.z_amino_acids./df_MAGs.Genome_size*1e6), 1, n_MAGs)
z_fattys        = reshape(convert(Array{Float64,1}, df_MAGs.z_fatty_acids./df_MAGs.Genome_size*1e6), 1, n_MAGs)
z_nucleos       = reshape(convert(Array{Float64,1}, df_MAGs.z_nucleotides./df_MAGs.Genome_size*1e6), 1, n_MAGs)
z_auxins        = reshape(convert(Array{Float64,1}, df_MAGs.z_auxins./df_MAGs.Genome_size*1e6), 1, n_MAGs)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

y_EM            = ones(size(V_cell,1))
########################################

########################################
# calc transporter density for glucose only
rho_ps          = zeros(1, size(V_cell,1))
y_DEs           = zeros(1, size(V_cell,1))
failed_MAGs     = Int[]

j               = glucose_idx
t_start         = time()

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
        y_DEs[1,i]  = 0.0
        push!(failed_MAGs, i)
    end

    # Progress update every 100 MAGs with time estimate
    if i % 100 == 0
        elapsed     = time() - t_start
        rate        = i / elapsed          # MAGs per second
        remaining   = (n_MAGs - i) / rate  # seconds remaining
        eta_min     = remaining / 60
        println("MAG $i of $n_MAGs | $(length(failed_MAGs)) failed | elapsed: $(round(elapsed/60, digits=1)) min | ETA: $(round(eta_min, digits=1)) min")
        # Intermediate save every 500 MAGs
        if i % 500 == 0
            save("files/my_MAGs_assimilation_glucose_partial.jld", "rho", rho_ps, "yDE", y_DEs, "progress", i)
        end
    end
end

println("\nCompleted all $n_MAGs MAGs")
println("Total convergence failures: $(length(failed_MAGs))")
if length(failed_MAGs) > 0
    println("Failed MAG indices: $failed_MAGs")
end

# Clean up problematic values
rho_ps[rho_ps.==0.0] .= 1e-12
rho_ps[isnan.(rho_ps)] .= 1e-12
rho_ps[isinf.(rho_ps)] .= 1e-12
rho_ps[rho_ps.<0.0]   .= 1e-12
y_DEs[isnan.(y_DEs)]  .= 0.0
y_DEs[isinf.(y_DEs)]  .= 0.0
println("Median transporter density: $(median(rho_ps))")

# N_C for glucose
N_C     = zeros(1)
N_C[1]  = DEBmicroTrait.extract_composition(df_metabolites.Formula[glucose_idx])[1]
println("N_C for glucose: $(N_C[1])")

N_SB    = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, rho_ps, Min_gen_time, Gram_stain)
Vmax    = @. 180.0*60^2*N_SB.*N_C
D_S     = DEBmicroTrait.aqueous_diffusivity([df_metabolites.Molecular_weight[glucose_idx]])
K_D     = DEBmicroTrait.specific_reference_affinity(V_cell, rho_ps, D_S)
a_s     = Vmax./K_D

########################################
# Final save
save("files/my_MAGs_assimilation_glucose.jld", "rho", rho_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C, "failed", failed_MAGs)
println("Saved: files/my_MAGs_assimilation_glucose.jld")
t_total = (time() - t_start) / 60
println("Total runtime: $(round(t_total, digits=1)) minutes")
########################################
