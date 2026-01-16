"""
    calc_attack.jl

Core optimization function for computing worst-case line flows or voltage magnitudes
given a set of compromised generators. Julia port of MATLAB calc_attack.m.

Uses Dan's exact power flow formulation with complex admittance arithmetic.
"""

using JuMP
using Ipopt
using PowerModels

include("attack_utils.jl")

"""
    calc_attack(data, wft, wtf, wv, compromised, alpha, pvpq_k;
                optimizer=Ipopt.Optimizer, verbose=false)

Compute worst-case line flows or voltage magnitudes for a given set of
compromised generators.

# Arguments
- `data`: PowerModels network data dictionary (solved OPF to get nominal setpoints)
- `wft`: Weight vector for maximizing I²_ft (from→to direction), length nbranch
- `wtf`: Weight vector for maximizing I²_tf (to→from direction), length nbranch
- `wv`: Weight vector for maximizing/minimizing voltage magnitudes, length nbus
- `compromised`: Vector of generator indices (1-based) that are compromised
- `alpha`: Vector of participation factors for AGC, length ngen
- `pvpq_k`: Slope parameter for smoothed PV/PQ switching characteristic

# Returns
Dict with optimization results.

# Notes
- Follows Dan's MATLAB calc_attack.m formulation exactly
- Assumes Pmax for compromised generators represents Smax (apparent power limit)
- Assumes input data is a solved OPF (provides nominal Pg, Qg, Vstar setpoints)
"""
function calc_attack(data::Dict{String,Any},
                     wft::Vector{Float64},
                     wtf::Vector{Float64},
                     wv::Vector{Float64},
                     compromised::Vector{Int},
                     alpha::Vector{Float64},
                     pvpq_k::Float64;
                     optimizer=Ipopt.Optimizer,
                     verbose::Bool=false)

    # Configuration (matching Dan's code)
    V_lv = 0.6  # Lower bound on voltage magnitudes
    pvpq_switching_on = true
    distslack_threshold_on = true

    # Use make_basic_network to get consecutive indexing (like MATPOWER ext2int)
    mpc = PowerModels.make_basic_network(deepcopy(data))

    # Extract dimensions
    nbus = length(mpc["bus"])
    ngen = length(mpc["gen"])
    nbranch = length(mpc["branch"])
    baseMVA = mpc["baseMVA"]

    # Get ordered keys (should now be 1:n after make_basic_network)
    bus_ids = sort([parse(Int, k) for k in keys(mpc["bus"])])
    gen_ids = sort([parse(Int, k) for k in keys(mpc["gen"])])
    branch_ids = sort([parse(Int, k) for k in keys(mpc["branch"])])

    # Create index mappings
    bus_id_to_idx = Dict(id => idx for (idx, id) in enumerate(bus_ids))

    # Extract bus data
    Pd = zeros(nbus)
    Qd = zeros(nbus)
    Gs = zeros(nbus)
    Bs = zeros(nbus)
    bus_type = zeros(Int, nbus)

    for (idx, bus_id) in enumerate(bus_ids)
        bus = mpc["bus"][string(bus_id)]
        # Shunts are on buses in PowerModels
        Gs[idx] = get(bus, "gs", 0.0) / baseMVA
        Bs[idx] = get(bus, "bs", 0.0) / baseMVA
        bus_type[idx] = get(bus, "bus_type", 1)
    end

    # CRITICAL: PowerModels stores loads in separate "load" dictionary, NOT on buses!
    # Aggregate all loads to their respective buses (like MATPOWER mpc.bus(:,PD))
    if haskey(mpc, "load")
        for (lk, load) in mpc["load"]
            load_bus = load["load_bus"]
            if haskey(bus_id_to_idx, load_bus)
                idx = bus_id_to_idx[load_bus]
                Pd[idx] += get(load, "pd", 0.0) / baseMVA
                Qd[idx] += get(load, "qd", 0.0) / baseMVA
            end
        end
    end

    # Find reference bus
    ref_bus_idx = findfirst(bus_type .== 3)
    if isnothing(ref_bus_idx)
        ref_bus_idx = 1
    end

    # Extract generator data
    Pg0 = zeros(ngen)
    Qg0 = zeros(ngen)
    Pgmax = zeros(ngen)
    Pgmin = zeros(ngen)
    Qgmax = zeros(ngen)
    Qgmin = zeros(ngen)
    Vstar = zeros(ngen)
    gen_bus_idx = zeros(Int, ngen)

    for (idx, gen_id) in enumerate(gen_ids)
        gen = mpc["gen"][string(gen_id)]
        Pg0[idx] = get(gen, "pg", 0.0) / baseMVA
        Qg0[idx] = get(gen, "qg", 0.0) / baseMVA
        Pgmax[idx] = get(gen, "pmax", 0.0) / baseMVA
        Pgmin[idx] = get(gen, "pmin", 0.0) / baseMVA
        Qgmax[idx] = get(gen, "qmax", 0.0) / baseMVA
        Qgmin[idx] = get(gen, "qmin", 0.0) / baseMVA
        Vstar[idx] = get(gen, "vg", 1.0)
        gen_bus = gen["gen_bus"]
        gen_bus_idx[idx] = bus_id_to_idx[gen_bus]
    end

    # Bump Q limits if Qg is near limit (Dan's code lines 81-82)
    for i in 1:ngen
        if Qgmin[i] != 0 && abs(Qg0[i] - Qgmin[i]) / abs(Qgmin[i]) <= 0.01
            Qgmin[i] *= 1.05
        end
        if Qgmax[i] != 0 && abs(Qg0[i] - Qgmax[i]) / abs(Qgmax[i]) <= 0.01
            Qgmax[i] *= 1.05
        end
    end

    # Create compromised boolean vector
    is_compromised = falses(ngen)
    for i in compromised
        if 1 <= i <= ngen
            is_compromised[i] = true
        end
    end

    # Adjust alpha for generators at limits (Dan's code line 87-88)
    alpha_adj = copy(alpha)
    for i in 1:ngen
        at_pmin = Pgmin[i] != 0 ? abs(Pg0[i] - Pgmin[i]) / abs(Pgmin[i]) <= 0.01 : false
        at_pmax = Pgmax[i] != 0 ? abs(Pg0[i] - Pgmax[i]) / abs(Pgmax[i]) <= 0.01 : false
        at_zero_limit = abs(Pg0[i] - Pgmin[i]) <= 0.01 && abs(Pgmin[i]) <= 0.01

        if at_pmin || at_pmax || at_zero_limit
            alpha_adj[i] = 0.0
        end
    end

    # Renormalize participation factors
    alpha_norm = sqrt(sum(alpha_adj.^2))
    if alpha_norm > 0
        alpha_adj ./= alpha_norm
    end

    # Extract branch data
    f_bus_idx = zeros(Int, nbranch)
    t_bus_idx = zeros(Int, nbranch)
    gl = zeros(nbranch)
    bl = zeros(nbranch)
    bsl = zeros(nbranch)
    tau = zeros(nbranch)
    theta_sh = zeros(nbranch)

    for (idx, branch_id) in enumerate(branch_ids)
        branch = mpc["branch"][string(branch_id)]
        f_bus_idx[idx] = bus_id_to_idx[branch["f_bus"]]
        t_bus_idx[idx] = bus_id_to_idx[branch["t_bus"]]

        # Line parameters (Dan's code lines 124-136)
        br_r = get(branch, "br_r", 0.0)
        br_x = get(branch, "br_x", 0.01)  # Avoid zero impedance
        zl = br_r + im * br_x
        yl = 1.0 / zl
        gl[idx] = real(yl)
        bl[idx] = imag(yl)
        bsl[idx] = get(branch, "br_b", 0.0)
        tau[idx] = get(branch, "tap", 1.0)
        if tau[idx] == 0
            tau[idx] = 1.0
        end
        theta_sh[idx] = get(branch, "shift", 0.0) * π / 180
    end

    # Create JuMP model
    model = Model(optimizer)
    if !verbose
        set_silent(model)
    end

    # Ipopt settings to match MATLAB fmincon
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "tol", 1e-5)
    set_optimizer_attribute(model, "constr_viol_tol", 1e-5)
    set_optimizer_attribute(model, "mu_strategy", "monotone")

    # Decision variables (Dan's code lines 94-98)
    @variable(model, V[i=1:nbus] >= V_lv, start = 1.0)  # No upper bound like Dan's code
    @variable(model, theta[1:nbus], start = 0.0)
    @variable(model, Pg[1:ngen])
    @variable(model, Qg[1:ngen])
    @variable(model, Delta)  # AGC deviation signal

    # Set starting points from nominal values
    for (idx, bus_id) in enumerate(bus_ids)
        bus = mpc["bus"][string(bus_id)]
        set_start_value(V[idx], get(bus, "vm", 1.0))
        set_start_value(theta[idx], get(bus, "va", 0.0))
    end
    for i in 1:ngen
        set_start_value(Pg[i], Pg0[i])
        set_start_value(Qg[i], Qg0[i])
    end
    set_start_value(Delta, 0.0)

    # Voltage at generator buses (Dan's code line 100)
    Vg = [V[gen_bus_idx[i]] for i in 1:ngen]

    # === Power Flow Equations (Dan's exact formulation, lines 103-146) ===
    # Using complex arithmetic notation from Dan's code:
    # Sft = (conj(yl)-1i*bsl/2)./tau.^2 .* Wf - conj(yl./tl) .* Wft
    # Stf = (conj(yl)-1i*bsl/2)         .* Wt - conj(yl)./tl .* Wtf

    # Branch flow expressions
    p_ft = Vector{Any}(undef, nbranch)
    q_ft = Vector{Any}(undef, nbranch)
    p_tf = Vector{Any}(undef, nbranch)
    q_tf = Vector{Any}(undef, nbranch)

    for k in 1:nbranch
        f = f_bus_idx[k]
        t = t_bus_idx[k]

        # Wf = V(f)^2, Wt = V(t)^2
        Wf = V[f]^2
        Wt = V[t]^2

        # Wft = V(f)*V(t)*exp(1i*(theta(f)-theta(t)))
        # Real part: V(f)*V(t)*cos(theta(f)-theta(t))
        # Imag part: V(f)*V(t)*sin(theta(f)-theta(t))
        theta_diff = theta[f] - theta[t]
        Wft_re = V[f] * V[t] * cos(theta_diff)
        Wft_im = V[f] * V[t] * sin(theta_diff)

        # tl = tau*exp(-1i*theta_sh), so 1/tl = (1/tau)*exp(1i*theta_sh)
        # conj(yl) = gl - 1i*bl
        # conj(yl/tl) = (gl - 1i*bl) * exp(-1i*theta_sh) / tau
        #             = (gl - 1i*bl) * (cos(theta_sh) - 1i*sin(theta_sh)) / tau
        # Real: (gl*cos(theta_sh) - bl*sin(theta_sh)) / tau
        # Imag: (-gl*sin(theta_sh) - bl*cos(theta_sh)) / tau

        cos_sh = cos(theta_sh[k])
        sin_sh = sin(theta_sh[k])

        # Sft = (gl - 1i*(bl + bsl/2)) / tau^2 * Wf - conj(yl/tl) * Wft
        # Real(Sft) = gl/tau^2 * Wf - Real(conj(yl/tl) * Wft)
        # Real(conj(yl/tl) * Wft) = Real((gl*cos_sh - bl*sin_sh)/tau - 1i*(gl*sin_sh + bl*cos_sh)/tau) * (Wft_re + 1i*Wft_im))
        #                        = (gl*cos_sh - bl*sin_sh)/tau * Wft_re + (gl*sin_sh + bl*cos_sh)/tau * Wft_im

        conj_yl_tl_re = (gl[k] * cos_sh - bl[k] * sin_sh) / tau[k]
        conj_yl_tl_im = -(gl[k] * sin_sh + bl[k] * cos_sh) / tau[k]

        # Real(conj_yl_tl * Wft) = conj_yl_tl_re * Wft_re - conj_yl_tl_im * Wft_im
        # Imag(conj_yl_tl * Wft) = conj_yl_tl_re * Wft_im + conj_yl_tl_im * Wft_re

        p_ft[k] = @expression(model,
            gl[k] / tau[k]^2 * Wf - (conj_yl_tl_re * Wft_re - conj_yl_tl_im * Wft_im)
        )

        q_ft[k] = @expression(model,
            -(bl[k] + bsl[k]/2) / tau[k]^2 * Wf - (conj_yl_tl_re * Wft_im + conj_yl_tl_im * Wft_re)
        )

        # Stf = (gl - 1i*(bl + bsl/2)) * Wt - conj(yl)/tl * Wtf
        # Wtf = conj(Wft) = Wft_re - 1i*Wft_im
        # conj(yl)/tl = (gl - 1i*bl) / (tau*exp(-1i*theta_sh))
        #             = (gl - 1i*bl) * exp(1i*theta_sh) / tau
        #             = (gl - 1i*bl) * (cos_sh + 1i*sin_sh) / tau
        # Real: (gl*cos_sh + bl*sin_sh) / tau
        # Imag: (gl*sin_sh - bl*cos_sh) / tau

        conj_yl_div_tl_re = (gl[k] * cos_sh + bl[k] * sin_sh) / tau[k]
        conj_yl_div_tl_im = (gl[k] * sin_sh - bl[k] * cos_sh) / tau[k]

        # Real(conj_yl_div_tl * Wtf) = conj_yl_div_tl_re * Wft_re + conj_yl_div_tl_im * Wft_im
        # Imag(conj_yl_div_tl * Wtf) = -conj_yl_div_tl_re * Wft_im + conj_yl_div_tl_im * Wft_re

        p_tf[k] = @expression(model,
            gl[k] * Wt - (conj_yl_div_tl_re * Wft_re + conj_yl_div_tl_im * Wft_im)
        )

        q_tf[k] = @expression(model,
            -(bl[k] + bsl[k]/2) * Wt - (-conj_yl_div_tl_re * Wft_im + conj_yl_div_tl_im * Wft_re)
        )
    end

    # === Power Balance Constraints (Dan's code lines 192-216) ===
    # P = sum(p_ft at bus) + sum(p_tf at bus) + V^2 * Gs
    # Q = sum(q_ft at bus) + sum(q_tf at bus) - V^2 * Bs
    # Pgsum - Pd == P
    # Qgsum - Qd == Q

    # Build lookup tables for which branches connect to each bus
    bus_arcs_from = [Int[] for _ in 1:nbus]
    bus_arcs_to = [Int[] for _ in 1:nbus]
    for k in 1:nbranch
        push!(bus_arcs_from[f_bus_idx[k]], k)
        push!(bus_arcs_to[t_bus_idx[k]], k)
    end

    # Build lookup for generators at each bus
    bus_gens = [Int[] for _ in 1:nbus]
    for i in 1:ngen
        push!(bus_gens[gen_bus_idx[i]], i)
    end

    # Power balance constraints
    for i in 1:nbus
        # Total P injection at bus
        P_inj = @expression(model,
            sum(p_ft[k] for k in bus_arcs_from[i]; init=0.0) +
            sum(p_tf[k] for k in bus_arcs_to[i]; init=0.0) +
            V[i]^2 * Gs[i]
        )

        # Total Q injection at bus
        Q_inj = @expression(model,
            sum(q_ft[k] for k in bus_arcs_from[i]; init=0.0) +
            sum(q_tf[k] for k in bus_arcs_to[i]; init=0.0) -
            V[i]^2 * Bs[i]
        )

        # Generator power at bus
        Pg_sum = @expression(model, sum(Pg[g] for g in bus_gens[i]; init=0.0))
        Qg_sum = @expression(model, sum(Qg[g] for g in bus_gens[i]; init=0.0))

        # Power balance
        @constraint(model, Pg_sum - Pd[i] == P_inj)
        @constraint(model, Qg_sum - Qd[i] == Q_inj)
    end

    # Angle reference constraint
    @constraint(model, theta[ref_bus_idx] == 0)

    # === Generator Constraints (Dan's code lines 220-315) ===

    for i in 1:ngen
        if is_compromised[i]
            # Compromised generators: box constraints + apparent power limit (lines 305-314)
            @constraint(model, Pg[i] >= Pgmin[i])
            @constraint(model, Pg[i] <= Pgmax[i])
            @constraint(model, Qg[i] >= Qgmin[i])
            @constraint(model, Qg[i] <= Qgmax[i])
            # Smax constraint (assumes Pmax = Smax for compromised gens)
            @constraint(model, Pg[i]^2 + Qg[i]^2 <= Pgmax[i]^2)
        else
            # Non-compromised generators

            # Active power: S-curve AGC model (lines 220-265)
            if distslack_threshold_on && alpha_adj[i] > 0
                # S-curve model for distributed slack
                k_scurve, Delta0 = optimal_scurve_k_asym_brute_force(
                    Pgmin[i], Pgmax[i], Pg0[i], alpha_adj[i])

                @constraint(model,
                    Pg[i] == Pgmin[i] + (Pgmax[i] - Pgmin[i]) /
                             (1 + exp(-k_scurve * (Delta - Delta0))))
            else
                # Fixed output (not participating in AGC)
                @constraint(model, Pg[i] == Pg0[i])
            end

            # Reactive power: S-curve PV/PQ switching model (lines 267-299)
            if pvpq_switching_on
                # Check if we can compute beta (avoid log of non-positive number)
                Qg_range = Qgmax[i] - Qgmin[i]
                numerator = Qgmax[i] - Qg0[i]
                denominator = Qg0[i] - Qgmin[i]

                if Qg_range > 1e-6 && numerator > 1e-6 && denominator > 1e-6
                    beta = log(numerator / denominator)

                    @constraint(model,
                        Qg[i] == Qgmin[i] + (Qgmax[i] - Qgmin[i]) /
                                 (1 + exp(pvpq_k * (Vg[i] - Vstar[i]) + beta)))
                else
                    # Degenerate case: fix Qg to nominal
                    @constraint(model, Qg[i] == Qg0[i])
                end
            else
                # Fixed voltage (PV bus model)
                @constraint(model, Vg[i] == Vstar[i])
            end
        end
    end

    # === Objective: Maximize Weighted Sum of I² and V (Dan's code lines 326-336) ===
    # Dan's exact squared current formulas

    Ift_sq = Vector{Any}(undef, nbranch)
    Itf_sq = Vector{Any}(undef, nbranch)

    for k in 1:nbranch
        f = f_bus_idx[k]
        t = t_bus_idx[k]
        theta_diff = theta[f] - theta[t]

        # Dan's formula (lines 326-332):
        # Ift_sq = V(f)^2/tau^4 * (gl^2 + (bl+bsl/2)^2)
        #        + V(t)^2/tau^2 * (gl^2 + bl^2)
        #        - 2*V(f)*V(t)/tau^3 * ((gl^2 + bl*(bl+bsl/2))*cos(θf-θt-θsh) - gl*bsl/2*sin(θf-θt-θsh))

        Ift_sq[k] = @expression(model,
            V[f]^2 / tau[k]^4 * (gl[k]^2 + (bl[k] + bsl[k]/2)^2) +
            V[t]^2 / tau[k]^2 * (gl[k]^2 + bl[k]^2) -
            2 * V[f] * V[t] / tau[k]^3 * (
                (gl[k]^2 + bl[k] * (bl[k] + bsl[k]/2)) * cos(theta_diff - theta_sh[k]) -
                gl[k] * bsl[k]/2 * sin(theta_diff - theta_sh[k])
            )
        )

        # Itf_sq = V(t)^2 * (gl^2 + (bl+bsl/2)^2)
        #        + V(f)^2/tau^2 * (gl^2 + bl^2)
        #        - 2*V(f)*V(t)/tau * ((gl^2 + bl*(bl+bsl/2))*cos(θf-θt-θsh) + gl*bsl/2*sin(θf-θt-θsh))

        Itf_sq[k] = @expression(model,
            V[t]^2 * (gl[k]^2 + (bl[k] + bsl[k]/2)^2) +
            V[f]^2 / tau[k]^2 * (gl[k]^2 + bl[k]^2) -
            2 * V[f] * V[t] / tau[k] * (
                (gl[k]^2 + bl[k] * (bl[k] + bsl[k]/2)) * cos(theta_diff - theta_sh[k]) +
                gl[k] * bsl[k]/2 * sin(theta_diff - theta_sh[k])
            )
        )
    end

    # Objective: maximize wft'*Ift_sq + wtf'*Itf_sq + wv'*V (Dan's line 336)
    # Note: Dan minimizes (-wft'*Ift_sq ...) which is equivalent to maximizing
    @objective(model, Max,
        sum(wft[k] * Ift_sq[k] for k in 1:nbranch) +
        sum(wtf[k] * Itf_sq[k] for k in 1:nbranch) +
        sum(wv[i] * V[i] for i in 1:nbus)
    )

    # Solve
    solve_start = time()
    optimize!(model)
    solve_time = time() - solve_start

    status = termination_status(model)
    success = status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL

    # Extract results
    result = Dict{String,Any}(
        "success" => success,
        "termination_status" => status,
        "solve_time" => solve_time,
        "objective" => success ? objective_value(model) : NaN,
        "V" => success ? value.(V) : fill(NaN, nbus),
        "theta" => success ? value.(theta) : fill(NaN, nbus),
        "Pg" => success ? value.(Pg) : fill(NaN, ngen),
        "Qg" => success ? value.(Qg) : fill(NaN, ngen),
        "Delta" => success ? value(Delta) : NaN,
        "Ift_sq" => success ? [value(Ift_sq[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "Itf_sq" => success ? [value(Itf_sq[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "p_ft" => success ? [value(p_ft[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "q_ft" => success ? [value(q_ft[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "p_tf" => success ? [value(p_tf[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "q_tf" => success ? [value(q_tf[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "nbus" => nbus,
        "ngen" => ngen,
        "nbranch" => nbranch,
        "baseMVA" => baseMVA
    )

    return result
end

"""
    calc_attack_line_ft(data, line_idx, compromised, alpha, pvpq_k; kwargs...)

Convenience function to maximize current flow on a specific line (from→to direction).
"""
function calc_attack_line_ft(data::Dict{String,Any}, line_idx::Int,
                             compromised::Vector{Int}, alpha::Vector{Float64},
                             pvpq_k::Float64; kwargs...)
    mpc = PowerModels.make_basic_network(deepcopy(data))
    nbranch = length(mpc["branch"])
    nbus = length(mpc["bus"])

    wft = zeros(nbranch)
    wtf = zeros(nbranch)
    wv = zeros(nbus)
    wft[line_idx] = 1.0
    return calc_attack(data, wft, wtf, wv, compromised, alpha, pvpq_k; kwargs...)
end

"""
    calc_attack_line_tf(data, line_idx, compromised, alpha, pvpq_k; kwargs...)

Convenience function to maximize current flow on a specific line (to→from direction).
"""
function calc_attack_line_tf(data::Dict{String,Any}, line_idx::Int,
                             compromised::Vector{Int}, alpha::Vector{Float64},
                             pvpq_k::Float64; kwargs...)
    mpc = PowerModels.make_basic_network(deepcopy(data))
    nbranch = length(mpc["branch"])
    nbus = length(mpc["bus"])

    wft = zeros(nbranch)
    wtf = zeros(nbranch)
    wv = zeros(nbus)
    wtf[line_idx] = 1.0
    return calc_attack(data, wft, wtf, wv, compromised, alpha, pvpq_k; kwargs...)
end

"""
    calc_attack_voltage_min(data, bus_idx, compromised, alpha, pvpq_k; kwargs...)

Convenience function to minimize voltage at a specific bus.
"""
function calc_attack_voltage_min(data::Dict{String,Any}, bus_idx::Int,
                                 compromised::Vector{Int}, alpha::Vector{Float64},
                                 pvpq_k::Float64; kwargs...)
    mpc = PowerModels.make_basic_network(deepcopy(data))
    nbranch = length(mpc["branch"])
    nbus = length(mpc["bus"])

    wft = zeros(nbranch)
    wtf = zeros(nbranch)
    wv = zeros(nbus)
    wv[bus_idx] = -1.0  # Negative weight to minimize
    return calc_attack(data, wft, wtf, wv, compromised, alpha, pvpq_k; kwargs...)
end

"""
    calc_attack_voltage_max(data, bus_idx, compromised, alpha, pvpq_k; kwargs...)

Convenience function to maximize voltage at a specific bus.
"""
function calc_attack_voltage_max(data::Dict{String,Any}, bus_idx::Int,
                                 compromised::Vector{Int}, alpha::Vector{Float64},
                                 pvpq_k::Float64; kwargs...)
    mpc = PowerModels.make_basic_network(deepcopy(data))
    nbranch = length(mpc["branch"])
    nbus = length(mpc["bus"])

    wft = zeros(nbranch)
    wtf = zeros(nbranch)
    wv = zeros(nbus)
    wv[bus_idx] = 1.0
    return calc_attack(data, wft, wtf, wv, compromised, alpha, pvpq_k; kwargs...)
end
