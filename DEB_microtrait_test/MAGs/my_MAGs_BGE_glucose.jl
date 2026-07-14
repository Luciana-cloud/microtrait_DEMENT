using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations

########################################
# I/O
df_MAGs           = CSV.read("files/MAGs2traits_DEBinput.csv", DataFrame, missingstring="N/A")
n_MAGs            = nrow(df_MAGs)
println("Loaded $n_MAGs MAGs")

# Load pre-processed parameter files
assimilation      = load("files/my_MAGs_assimilation_glucose.jld")
enzymes           = load("files/my_MAGs_enzymes.jld")
maintenance       = load("files/my_MAGs_maintenance.jld")
protein_synthesis = load("files/my_MAGs_protein_synthesis.jld")
turnover          = load("files/my_MAGs_turnover.jld")
initb             = load("files/my_MAGs_batch_init.jld")

# Glucose index in our glucose-only assimilation file
glucose_idx = 1
println("All parameter files loaded successfully")
########################################

########################################
# Callback to stop simulation when substrate is depleted
condition(u,t,integrator) = u[1] - 1e-5
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition, affect!)

# Pre-allocate time series arrays
BGE_tseries     = zeros(n_MAGs, 500)
BR_tseries      = zeros(n_MAGs, 500)
BP_tseries      = zeros(n_MAGs, 500)
r_tseries       = zeros(n_MAGs, 500)
failed_sims     = Int[]

# Pre-allocate median arrays
BGE_median      = fill(NaN, n_MAGs)
BP_median       = fill(NaN, n_MAGs)
BR_median       = fill(NaN, n_MAGs)
r_median        = fill(NaN, n_MAGs)

t_start         = time()
########################################

########################################
# Helper function to save partial CSV
function save_partial_csv(up_to::Int)
    df_partial              = DataFrame()
    df_partial.MAG_id       = df_MAGs.Isolate[1:up_to]
    df_partial.genomesize   = df_MAGs.Genome_size[1:up_to]
    df_partial.mingt        = df_MAGs.Min_gen_time[1:up_to]
    df_partial.rrn          = df_MAGs.rRNA_genes[1:up_to]
    df_partial.BGE          = BGE_median[1:up_to]
    df_partial.rgrowth      = r_median[1:up_to]
    df_partial.BP           = BP_median[1:up_to]
    df_partial.BR           = BR_median[1:up_to]
    CSV.write("files/my_MAGs_BGE_glucose_partial.csv", df_partial)
end
########################################

########################################
# Run simulation for each MAG
for i in 1:n_MAGs
    try
        id_MAG     = i
        p          = DEBmicroTrait.init_batch_model(id_MAG, glucose_idx, assimilation, enzymes, maintenance, protein_synthesis, turnover)
        n_polymers = p.setup_pars.n_polymers
        n_monomers = p.setup_pars.n_monomers
        n_microbes = p.setup_pars.n_microbes

        u0 = zeros(p.setup_pars.dim)
        u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_MAG]
        u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_MAG]
        u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.25

        tspan = (0.0, 1000.0)
        prob  = ODEProblem(DEBmicroTrait.batch_model!, u0, tspan, p)
        sol   = solve(prob, alg_hints=[:stiff], callback=cb)

        du  = zeros(p.setup_pars.dim)
        BR  = [DEBmicroTrait.batch_model!(du, sol.u[k], p, 0)[end] for k in 1:size(sol.t,1)]
        BP  = [DEBmicroTrait.batch_model!(du, sol.u[k], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[k], p, 0)[3] for k in 1:size(sol.t,1)]
        BGE = @. BP/(BP + BR)
        r   = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[k][2]], [sol[k][3]])[1] for k in 1:size(sol.t,1)]

        for k in 1:length(sol.t)
            BGE_tseries[i,k] = BGE[k]
            BR_tseries[i,k]  = BR[k]
            BP_tseries[i,k]  = BP[k]
            r_tseries[i,k]   = r[k]
        end

        # Compute median for this MAG immediately
        BGE_clean = filter(x -> x > 0.0 && x <= 1.0, BGE_tseries[i,:])
        BGE_median[i] = length(BGE_clean) > 0 ? median(BGE_clean) : NaN

        BP_clean = filter(!iszero, BP_tseries[i,:])
        BP_median[i] = length(BP_clean) > 0 ? median(BP_clean) : NaN

        BR_clean = filter(!iszero, BR_tseries[i,:])
        BR_median[i] = length(BR_clean) > 0 ? median(BR_clean) : NaN

        r_clean = filter(!iszero, r_tseries[i,:])
        r_median[i] = length(r_clean) > 0 ? median(r_clean) : NaN

    catch e
        push!(failed_sims, i)
    end

    # Progress update every 100 MAGs with ETA
    if i % 100 == 0
        elapsed   = time() - t_start
        rate      = i / elapsed
        remaining = (n_MAGs - i) / rate
        eta_min   = remaining / 60
        println("MAG $i of $n_MAGs | $(length(failed_sims)) failed | elapsed: $(round(elapsed/60, digits=1)) min | ETA: $(round(eta_min, digits=1)) min")
    end

    # Save partial CSV every 500 MAGs
    if i % 500 == 0
        save_partial_csv(i)
        println("  Partial CSV saved up to MAG $i")
    end
end

println("\nSimulation complete!")
println("Total failed simulations: $(length(failed_sims))")
if length(failed_sims) > 0
    println("Failed MAG indices: $failed_sims")
end
########################################

########################################
# Build final output DataFrame
df_out              = DataFrame()
df_out.MAG_id       = df_MAGs.Isolate
df_out.genomesize   = df_MAGs.Genome_size
df_out.mingt        = df_MAGs.Min_gen_time
df_out.rrn          = df_MAGs.rRNA_genes
df_out.BGE          = BGE_median
df_out.rgrowth      = r_median
df_out.BP           = BP_median
df_out.BR           = BR_median

# Save final CSV
CSV.write("files/my_MAGs_BGE_glucose.csv", df_out)
println("Saved: files/my_MAGs_BGE_glucose.csv")

# Summary statistics
valid_BGE = filter(!isnan, BGE_median)
println("\nSummary:")
println("MAGs with valid BGE: $(length(valid_BGE)) of $n_MAGs")
println("Median simulated BGE (CUE): $(median(valid_BGE))")
println("Min BGE: $(minimum(valid_BGE))")
println("Max BGE: $(maximum(valid_BGE))")
t_total = (time() - t_start) / 60
println("Total runtime: $(round(t_total, digits=1)) minutes")
########################################
