"""
    run_calc_attack.jl

High-level orchestration script for computing worst-case attack outcomes.
Julia port of MATLAB run_calc_attack.m.

Iterates over all lines (both flow directions) and all buses (min/max voltage)
to find the worst-case impacts of compromised generators.

Usage:
    julia --project=.. run_calc_attack.jl

Or from Julia REPL:
    include("run_calc_attack.jl")
    results = run_attack_analysis("../hawaii40.m", collect(ngen-9:ngen))
"""

using PowerModels
using JuMP
using Ipopt
using Printf

include("calc_attack.jl")

"""
    run_attack_analysis(case_path, compromised_gens;
                        pvpq_k=100.0, pf=0.95, optimizer=Ipopt.Optimizer, verbose=false)

Run complete attack analysis on a power system case.

# Arguments
- `case_path`: Path to MATPOWER case file
- `compromised_gens`: Vector of generator indices that are compromised

# Keyword Arguments
- `pvpq_k`: Slope for PV/PQ switching S-curve (default: 100.0)
- `pf`: Power factor for load Qd auto-fix (default: 0.95)
- `optimizer`: JuMP optimizer (default: Ipopt.Optimizer)
- `verbose`: Show solver output (default: false)

# Returns
Dict with:
- `Ift`: Worst-case current magnitudes (from→to), length nbranch
- `Itf`: Worst-case current magnitudes (to→from), length nbranch
- `Vlow`: Worst-case minimum voltages, length nbus
- `Vhigh`: Worst-case maximum voltages, length nbus
- `baseline`: Baseline OPF results
- `case_path`: Input case path
- `compromised_gens`: Compromised generator indices
"""
function run_attack_analysis(case_path::String, compromised_gens::Vector{Int};
                             pvpq_k::Float64=100.0, pf::Float64=0.95,
                             optimizer=Ipopt.Optimizer, verbose::Bool=false)

    println("=" ^ 60)
    println("Attack Analysis: ", case_path)
    println("=" ^ 60)

    # Load case
    println("Loading case...")
    mpc = PowerModels.parse_file(case_path)

    # Remove generators that are turned off
    gens_to_remove = String[]
    for (gk, gen) in mpc["gen"]
        if get(gen, "gen_status", 1) == 0
            push!(gens_to_remove, gk)
        end
    end
    for gk in gens_to_remove
        delete!(mpc["gen"], gk)
    end
    println("Removed $(length(gens_to_remove)) offline generators")

    # Fix Qd if all zero (before building ref)
    all_qd_zero = true
    for (lk, load) in mpc["load"]
        if get(load, "qd", 0.0) != 0.0
            all_qd_zero = false
            break
        end
    end
    if all_qd_zero
        println("Auto-fixing Qd with power factor $(pf)")
        for (lk, load) in mpc["load"]
            pd = get(load, "pd", 0.0)
            load["qd"] = pd * sqrt((1/pf^2) - 1)
        end
    end

    # Get dimensions using build_ref
    pm_ref = PowerModels.build_ref(mpc)[:it][:pm][:nw][0]
    nbus = length(pm_ref[:bus])
    ngen = length(pm_ref[:gen])
    nbranch = length(pm_ref[:branch])
    baseMVA = pm_ref[:baseMVA]

    println("Network: $(nbus) buses, $(ngen) generators, $(nbranch) branches")
    println("Compromised generators: ", compromised_gens)

    # Set lower generation bounds of compromised generators to zero
    # (models ability to shut off generators)
    gen_ids = sort(collect(keys(pm_ref[:gen])))
    for (idx, gen_id) in enumerate(gen_ids)
        if idx in compromised_gens
            mpc["gen"][string(gen_id)]["pmin"] = 0.0
        end
    end

    # Run baseline OPF to get nominal values
    println("Running baseline OPF...")
    baseline = PowerModels.solve_opf(mpc, ACPPowerModel, optimizer)

    if baseline["termination_status"] != MOI.LOCALLY_SOLVED &&
       baseline["termination_status"] != MOI.OPTIMAL
        error("Baseline OPF failed: $(baseline["termination_status"])")
    end
    println("Baseline OPF solved successfully")

    # Update mpc with baseline solution
    for (bk, bus_sol) in baseline["solution"]["bus"]
        if haskey(mpc["bus"], bk)
            mpc["bus"][bk]["vm"] = bus_sol["vm"]
            mpc["bus"][bk]["va"] = bus_sol["va"]
        end
    end
    for (gk, gen_sol) in baseline["solution"]["gen"]
        if haskey(mpc["gen"], gk)
            mpc["gen"][gk]["pg"] = gen_sol["pg"]
            mpc["gen"][gk]["qg"] = gen_sol["qg"]
            # Set vg from bus voltage
            gen_bus = mpc["gen"][gk]["gen_bus"]
            if haskey(baseline["solution"]["bus"], string(gen_bus))
                mpc["gen"][gk]["vg"] = baseline["solution"]["bus"][string(gen_bus)]["vm"]
            end
        end
    end

    # Rebuild ref after updating with baseline solution
    pm_ref = PowerModels.build_ref(mpc)[:it][:pm][:nw][0]
    gen_ids = sort(collect(keys(pm_ref[:gen])))

    # Compute participation factors based on nominal outputs
    Pg0 = zeros(ngen)
    for (idx, gen_id) in enumerate(gen_ids)
        gen = pm_ref[:gen][gen_id]
        Pg0[idx] = gen["pg"] / baseMVA
    end
    Pg_total = sum(Pg0)
    alpha = Pg_total > 0 ? Pg0 ./ Pg_total : zeros(ngen)

    # Preallocate results
    Ift = fill(NaN, nbranch)
    Itf = fill(NaN, nbranch)
    Vlow = fill(NaN, nbus)
    Vhigh = fill(NaN, nbus)

    # Results storage for detailed analysis
    results_ft = Vector{Dict{String,Any}}(undef, nbranch)
    results_tf = Vector{Dict{String,Any}}(undef, nbranch)
    results_vlow = Vector{Dict{String,Any}}(undef, nbus)
    results_vhigh = Vector{Dict{String,Any}}(undef, nbus)

    println("\n--- Computing Worst-Case Line Flows (from→to) ---")
    for i in 1:nbranch
        @printf("Line %d/%d (from→to)...", i, nbranch)
        res = calc_attack_line_ft(mpc, i, compromised_gens, alpha, pvpq_k;
                                  optimizer=optimizer, verbose=verbose)
        results_ft[i] = res
        if res["success"]
            Ift[i] = sqrt(max(res["objective"], 0.0))
            @printf(" I = %.4f p.u.\n", Ift[i])
        else
            @printf(" FAILED (%s)\n", res["termination_status"])
        end
    end

    println("\n--- Computing Worst-Case Line Flows (to→from) ---")
    for i in 1:nbranch
        @printf("Line %d/%d (to→from)...", i, nbranch)
        res = calc_attack_line_tf(mpc, i, compromised_gens, alpha, pvpq_k;
                                  optimizer=optimizer, verbose=verbose)
        results_tf[i] = res
        if res["success"]
            Itf[i] = sqrt(max(res["objective"], 0.0))
            @printf(" I = %.4f p.u.\n", Itf[i])
        else
            @printf(" FAILED (%s)\n", res["termination_status"])
        end
    end

    println("\n--- Computing Worst-Case Minimum Voltages ---")
    for i in 1:nbus
        @printf("Bus %d/%d (min V)...", i, nbus)
        res = calc_attack_voltage_min(mpc, i, compromised_gens, alpha, pvpq_k;
                                      optimizer=optimizer, verbose=verbose)
        results_vlow[i] = res
        if res["success"]
            # Objective is negative of V for minimization
            Vlow[i] = -res["objective"]
            @printf(" V = %.4f p.u.\n", Vlow[i])
        else
            @printf(" FAILED (%s)\n", res["termination_status"])
        end
    end

    println("\n--- Computing Worst-Case Maximum Voltages ---")
    for i in 1:nbus
        @printf("Bus %d/%d (max V)...", i, nbus)
        res = calc_attack_voltage_max(mpc, i, compromised_gens, alpha, pvpq_k;
                                      optimizer=optimizer, verbose=verbose)
        results_vhigh[i] = res
        if res["success"]
            Vhigh[i] = res["objective"]
            @printf(" V = %.4f p.u.\n", Vhigh[i])
        else
            @printf(" FAILED (%s)\n", res["termination_status"])
        end
    end

    # Summary
    println("\n" * "=" ^ 60)
    println("SUMMARY")
    println("=" ^ 60)

    valid_ift = filter(!isnan, Ift)
    valid_itf = filter(!isnan, Itf)
    valid_vlow = filter(!isnan, Vlow)
    valid_vhigh = filter(!isnan, Vhigh)

    if !isempty(valid_ift)
        @printf("Max line current (f→t): %.4f p.u. at line %d\n",
                maximum(valid_ift), argmax(Ift))
    end
    if !isempty(valid_itf)
        @printf("Max line current (t→f): %.4f p.u. at line %d\n",
                maximum(valid_itf), argmax(Itf))
    end
    if !isempty(valid_vlow)
        @printf("Min voltage: %.4f p.u. at bus %d\n",
                minimum(valid_vlow), argmin(Vlow))
    end
    if !isempty(valid_vhigh)
        @printf("Max voltage: %.4f p.u. at bus %d\n",
                maximum(valid_vhigh), argmax(Vhigh))
    end

    @printf("\nSuccess rate: %d/%d (%.1f%%)\n",
            length(valid_ift) + length(valid_itf) + length(valid_vlow) + length(valid_vhigh),
            2*nbranch + 2*nbus,
            100 * (length(valid_ift) + length(valid_itf) + length(valid_vlow) + length(valid_vhigh)) / (2*nbranch + 2*nbus))

    return Dict{String,Any}(
        "Ift" => Ift,
        "Itf" => Itf,
        "Vlow" => Vlow,
        "Vhigh" => Vhigh,
        "baseline" => baseline,
        "case_path" => case_path,
        "compromised_gens" => compromised_gens,
        "nbus" => nbus,
        "ngen" => ngen,
        "nbranch" => nbranch,
        "results_ft" => results_ft,
        "results_tf" => results_tf,
        "results_vlow" => results_vlow,
        "results_vhigh" => results_vhigh
    )
end

"""
    print_voltage_violations(results; v_min=0.9, v_max=1.1)

Print summary of buses with voltage violations under worst-case attack.
"""
function print_voltage_violations(results::Dict{String,Any};
                                  v_min::Float64=0.9, v_max::Float64=1.1)
    Vlow = results["Vlow"]
    Vhigh = results["Vhigh"]

    println("\n--- Voltage Violations ---")
    println("Buses with worst-case low voltage < $(v_min) p.u.:")
    for i in 1:length(Vlow)
        if !isnan(Vlow[i]) && Vlow[i] < v_min
            @printf("  Bus %d: %.4f p.u.\n", i, Vlow[i])
        end
    end

    println("\nBuses with worst-case high voltage > $(v_max) p.u.:")
    for i in 1:length(Vhigh)
        if !isnan(Vhigh[i]) && Vhigh[i] > v_max
            @printf("  Bus %d: %.4f p.u.\n", i, Vhigh[i])
        end
    end
end

# Main script execution when run directly
if abspath(PROGRAM_FILE) == @__FILE__
    # Default: run on Hawaii40 case with last 10 generators compromised
    case_path = joinpath(@__DIR__, "..", "hawaii40.m")

    if !isfile(case_path)
        error("Case file not found: $case_path")
    end

    # Load case to determine number of generators
    mpc_temp = PowerModels.parse_file(case_path)
    ngen_temp = length(mpc_temp["gen"])

    # Compromise last 10 generators (or all if fewer than 10)
    n_compromised = min(10, ngen_temp)
    compromised_gens = collect((ngen_temp - n_compromised + 1):ngen_temp)

    println("Running attack analysis on Hawaii40 case")
    println("Compromising generators: ", compromised_gens)

    results = run_attack_analysis(case_path, compromised_gens; verbose=false)

    print_voltage_violations(results)

    println("\nAnalysis complete!")
end
