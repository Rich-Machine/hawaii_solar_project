"""
    calc_attack.jl

Core optimization function for computing worst-case line flows or voltage magnitudes
given a set of compromised generators. Direct Julia port of MATLAB calc_attack.m.

Usage:
    include("calc_attack.jl")
    res = calc_attack(data, wft, wtf, wv, c, alpha, pvpq_k)
"""

using JuMP
using Ipopt
using PowerModels

include("attack_utils.jl")

"""
    calc_attack(data, wft, wtf, wv, c, alpha, pvpq_k; optimizer, verbose)

Compute worst-case line flows or voltage magnitudes for a given set of compromised
generators. Direct port of Dan's MATLAB calc_attack.m.

Arguments:
- data:    PowerModels network data dictionary (with OPF solution as nominal setpoints)
- wft:     Weight vector for from->to current squared objective (length nbranch)
- wtf:     Weight vector for to->from current squared objective (length nbranch)
- wv:      Weight vector for voltage magnitudes objective (length nbus)
- c:       Boolean vector, true = compromised generator (length ngen)
- alpha:   Generator participation factors for AGC (length ngen)
- pvpq_k:  Slope parameter for smoothed PV/PQ switching S-curve

Returns:
Dict with fields:
- "success":  true if solver converged
- "f":        objective value (matching MATLAB res.f = double(-obj))
- "V":        voltage magnitudes (p.u.)
- "theta":    voltage angles (radians)
- "Pg":       active power outputs (p.u.)
- "Qg":       reactive power outputs (p.u.)
- "Delta":    AGC deviation signal
- "p_ft", "q_ft", "p_tf", "q_tf": branch power flows (p.u.)
"""
function calc_attack(data::Dict{String,Any},
                     wft::Vector{Float64},
                     wtf::Vector{Float64},
                     wv::Vector{Float64},
                     c::AbstractVector{Bool},
                     alpha::Vector{Float64},
                     pvpq_k::Float64;
                     optimizer=Ipopt.Optimizer,
                     verbose::Bool=false)

    # -------------------------------------------------------------------------
    # Configuration (matching Dan's MATLAB calc_attack.m)
    # -------------------------------------------------------------------------
    V_lv = 0.6           # Lower bound on voltage magnitudes
    pvpq_switching_on     = true
    distslack_threshold_on = true

    # Use make_basic_network to get consecutive indexing (equivalent to ext2int)
    mpc = PowerModels.make_basic_network(deepcopy(data))

    nbus    = length(mpc["bus"])
    ngen    = length(mpc["gen"])
    nbranch = length(mpc["branch"])

    bus_ids    = sort([parse(Int, k) for k in keys(mpc["bus"])])
    gen_ids    = sort([parse(Int, k) for k in keys(mpc["gen"])])
    branch_ids = sort([parse(Int, k) for k in keys(mpc["branch"])])

    bus_id_to_idx = Dict(id => idx for (idx, id) in enumerate(bus_ids))

    # -------------------------------------------------------------------------
    # Load mpc data (matching Dan's MATLAB lines 60-88)
    # -------------------------------------------------------------------------

    # Bus data
    bus_type = zeros(Int, nbus)
    Gs       = zeros(nbus)
    Bs       = zeros(nbus)
    for (idx, bus_id) in enumerate(bus_ids)
        bus = mpc["bus"][string(bus_id)]
        bus_type[idx] = get(bus, "bus_type", 1)
        Gs[idx] = get(bus, "gs", 0.0)
        Bs[idx] = get(bus, "bs", 0.0)
    end

    # Loads aggregated to buses (PowerModels stores loads separately from buses)
    Pd = zeros(nbus)
    Qd = zeros(nbus)
    if haskey(mpc, "load")
        for (_, load) in mpc["load"]
            load_bus = load["load_bus"]
            if haskey(bus_id_to_idx, load_bus)
                idx = bus_id_to_idx[load_bus]
                Pd[idx] += get(load, "pd", 0.0)
                Qd[idx] += get(load, "qd", 0.0)
            end
        end
    end

    # Reference bus (bus_type == 3)
    ref = findfirst(bus_type .== 3)
    if isnothing(ref)
        ref = 1
    end

    # Generator data
    Pg0     = zeros(ngen); Qg0    = zeros(ngen)
    Pgmax   = zeros(ngen); Pgmin  = zeros(ngen)
    Qgmax   = zeros(ngen); Qgmin  = zeros(ngen)
    Vstar   = zeros(ngen)
    gen_bus_idx = zeros(Int, ngen)

    for (idx, gen_id) in enumerate(gen_ids)
        gen = mpc["gen"][string(gen_id)]
        Pg0[idx]        = get(gen, "pg",   0.0)
        Qg0[idx]        = get(gen, "qg",   0.0)
        Pgmax[idx]      = get(gen, "pmax", 0.0)
        Pgmin[idx]      = get(gen, "pmin", 0.0)
        Qgmax[idx]      = get(gen, "qmax", 0.0)
        Qgmin[idx]      = get(gen, "qmin", 0.0)
        Vstar[idx]      = get(gen, "vg",   1.0)
        gen_bus_idx[idx] = bus_id_to_idx[gen["gen_bus"]]
    end

    # Bump Q limits when Qg is near its limit (Dan's MATLAB lines 81-82)
    for i in 1:ngen
        if Qgmin[i] != 0 && abs(Qg0[i] - Qgmin[i]) / abs(Qgmin[i]) <= 1e-2
            Qgmin[i] *= 1.05
        end
        if Qgmax[i] != 0 && abs(Qg0[i] - Qgmax[i]) / abs(Qgmax[i]) <= 1e-2
            Qgmax[i] *= 1.05
        end
    end

    # Zero out alpha for generators near their active power limits (Dan's MATLAB line 87)
    alpha_adj = copy(alpha)
    for i in 1:ngen
        at_pmin     = (Pgmin[i] != 0) && (abs(Pg0[i] - Pgmin[i]) / abs(Pgmin[i]) <= 1e-2)
        at_pmax     = (Pgmax[i] != 0) && (abs(Pg0[i] - Pgmax[i]) / abs(Pgmax[i]) <= 1e-2)
        at_zero_lim = (abs(Pg0[i] - Pgmin[i]) <= 1e-2) && (abs(Pgmin[i]) <= 1e-2)
        if at_pmin || at_pmax || at_zero_lim
            alpha_adj[i] = 0.0
        end
    end
    # Renormalize (Dan's MATLAB line 88)
    alpha_norm = sqrt(sum(alpha_adj.^2))
    if alpha_norm > 0
        alpha_adj ./= alpha_norm
    end

    # Branch data (Dan's MATLAB lines 106-141)
    f_bus_idx = zeros(Int, nbranch)
    t_bus_idx = zeros(Int, nbranch)
    gl        = zeros(nbranch)
    bl        = zeros(nbranch)
    bsl       = zeros(nbranch)
    tau_br    = zeros(nbranch)
    theta_sh  = zeros(nbranch)

    for (idx, branch_id) in enumerate(branch_ids)
        branch = mpc["branch"][string(branch_id)]
        f_bus_idx[idx] = bus_id_to_idx[branch["f_bus"]]
        t_bus_idx[idx] = bus_id_to_idx[branch["t_bus"]]
        br_r = get(branch, "br_r", 0.0)
        br_x = get(branch, "br_x", 0.01)
        zl   = br_r + im * br_x
        yl   = 1.0 / zl
        gl[idx]       = real(yl)
        bl[idx]       = imag(yl)
        bsl[idx]      = get(branch, "br_b",  0.0)
        tau_br[idx]   = get(branch, "tap",   1.0)
        if tau_br[idx] == 0; tau_br[idx] = 1.0; end
        theta_sh[idx] = get(branch, "shift", 0.0) * π / 180
    end

    # -------------------------------------------------------------------------
    # Create variables (Dan's MATLAB lines 94-100)
    # -------------------------------------------------------------------------
    model = Model(optimizer)
    if !verbose
        set_silent(model)
    end

    # Ipopt settings (approximating Dan's fmincon settings)
    set_optimizer_attribute(model, "max_iter",        1000)
    set_optimizer_attribute(model, "tol",             1e-5)
    set_optimizer_attribute(model, "constr_viol_tol", 1e-5)
    set_optimizer_attribute(model, "mu_strategy",     "monotone")

    @variable(model, V[i=1:nbus] >= V_lv, start = 1.0)
    @variable(model, theta[1:nbus],        start = 0.0)
    @variable(model, Pg[1:ngen])
    @variable(model, Qg[1:ngen])
    @variable(model, Delta)

    # Warm start from nominal OPF values
    for (idx, bus_id) in enumerate(bus_ids)
        bus = mpc["bus"][string(bus_id)]
        set_start_value(V[idx],     get(bus, "vm", 1.0))
        set_start_value(theta[idx], get(bus, "va", 0.0))
    end
    for i in 1:ngen
        set_start_value(Pg[i], Pg0[i])
        set_start_value(Qg[i], Qg0[i])
    end
    set_start_value(Delta, 0.0)

    Vg = [V[gen_bus_idx[i]] for i in 1:ngen]

    # -------------------------------------------------------------------------
    # Power flow expressions (Dan's MATLAB lines 103-193)
    # -------------------------------------------------------------------------
    p_ft = Vector{Any}(undef, nbranch)
    q_ft = Vector{Any}(undef, nbranch)
    p_tf = Vector{Any}(undef, nbranch)
    q_tf = Vector{Any}(undef, nbranch)

    for k in 1:nbranch
        f = f_bus_idx[k]
        t = t_bus_idx[k]
        Wf = V[f]^2
        Wt = V[t]^2
        θd = theta[f] - theta[t]
        Wft_re = V[f] * V[t] * cos(θd)
        Wft_im = V[f] * V[t] * sin(θd)

        cos_sh = cos(theta_sh[k])
        sin_sh = sin(theta_sh[k])
        τ      = tau_br[k]

        # conj(yl/tl) where tl = tau*exp(-i*theta_sh)
        cyl_tl_re =  (gl[k]*cos_sh - bl[k]*sin_sh) / τ
        cyl_tl_im = -(gl[k]*sin_sh + bl[k]*cos_sh) / τ

        # conj(yl)/tl
        cyl_div_tl_re = (gl[k]*cos_sh + bl[k]*sin_sh) / τ
        cyl_div_tl_im = (gl[k]*sin_sh - bl[k]*cos_sh) / τ

        # Sft = (conj(yl) - i*bsl/2)/tau^2 * Wf - conj(yl/tl) * Wft
        p_ft[k] = @expression(model,
            gl[k]/τ^2 * Wf - (cyl_tl_re*Wft_re - cyl_tl_im*Wft_im))
        q_ft[k] = @expression(model,
            -(bl[k]+bsl[k]/2)/τ^2 * Wf - (cyl_tl_re*Wft_im + cyl_tl_im*Wft_re))

        # Stf = (conj(yl) - i*bsl/2) * Wt - conj(yl)/tl * Wtf
        p_tf[k] = @expression(model,
            gl[k]*Wt - (cyl_div_tl_re*Wft_re + cyl_div_tl_im*Wft_im))
        q_tf[k] = @expression(model,
            -(bl[k]+bsl[k]/2)*Wt - (-cyl_div_tl_re*Wft_im + cyl_div_tl_im*Wft_re))
    end

    # -------------------------------------------------------------------------
    # Power balance constraints (Dan's MATLAB lines 150-216)
    # -------------------------------------------------------------------------
    bus_arcs_from = [Int[] for _ in 1:nbus]
    bus_arcs_to   = [Int[] for _ in 1:nbus]
    for k in 1:nbranch
        push!(bus_arcs_from[f_bus_idx[k]], k)
        push!(bus_arcs_to[t_bus_idx[k]],   k)
    end
    bus_gens = [Int[] for _ in 1:nbus]
    for i in 1:ngen
        push!(bus_gens[gen_bus_idx[i]], i)
    end

    for i in 1:nbus
        P_inj  = @expression(model,
            sum(p_ft[k] for k in bus_arcs_from[i]; init=0.0) +
            sum(p_tf[k] for k in bus_arcs_to[i];   init=0.0) +
            V[i]^2 * Gs[i])
        Q_inj  = @expression(model,
            sum(q_ft[k] for k in bus_arcs_from[i]; init=0.0) +
            sum(q_tf[k] for k in bus_arcs_to[i];   init=0.0) -
            V[i]^2 * Bs[i])
        Pg_sum = @expression(model, sum(Pg[g] for g in bus_gens[i]; init=0.0))
        Qg_sum = @expression(model, sum(Qg[g] for g in bus_gens[i]; init=0.0))
        @constraint(model, Pg_sum - Pd[i] == P_inj)
        @constraint(model, Qg_sum - Qd[i] == Q_inj)
    end

    # Angle reference (Dan's MATLAB line 201-203)
    @constraint(model, theta[ref] == 0)

    # -------------------------------------------------------------------------
    # Generator constraints (Dan's MATLAB lines 220-315)
    # -------------------------------------------------------------------------
    for i in 1:ngen
        if c[i]
            # Compromised generators: apparent power limit + Pmin only.
            # NOTE: Dan's MATLAB does NOT add Pmax, Qmin, Qmax box constraints
            # for compromised generators — only Smax and Pmin.
            @constraint(model, Pg[i]^2 + Qg[i]^2 <= Pgmax[i]^2)  # Smax
            @constraint(model, Pg[i] >= Pgmin[i])
        else
            # Non-compromised: S-curve AGC model for active power (lines 220-265)
            if distslack_threshold_on && alpha_adj[i] > 0
                k_sc, Delta0 = optimal_scurve_k_asym_brute_force(
                    Pgmin[i], Pgmax[i], Pg0[i], alpha_adj[i])
                @constraint(model,
                    Pg[i] == Pgmin[i] + (Pgmax[i] - Pgmin[i]) /
                             (1 + exp(-k_sc * (Delta - Delta0))))
            else
                # Generator not participating in AGC: fix to nominal output
                @constraint(model, Pg[i] == Pg0[i])
            end

            # Non-compromised: PV/PQ switching S-curve for reactive power (lines 267-299)
            if pvpq_switching_on
                Qrange = Qgmax[i] - Qgmin[i]
                num    = Qgmax[i] - Qg0[i]
                den    = Qg0[i]   - Qgmin[i]
                if Qrange > 1e-6 && num > 1e-6 && den > 1e-6
                    beta = log(num / den)
                    @constraint(model,
                        Qg[i] == Qgmin[i] + (Qgmax[i] - Qgmin[i]) /
                                 (1 + exp(pvpq_k * (Vg[i] - Vstar[i]) + beta)))
                else
                    @constraint(model, Qg[i] == Qg0[i])
                end
            else
                # PV bus model: fixed voltage
                @constraint(model, Vg[i] == Vstar[i])
            end
        end
    end

    # Lower voltage limit (line 319-321)
    # (already enforced via @variable lower bound V >= V_lv above)

    # -------------------------------------------------------------------------
    # Objective function (Dan's MATLAB lines 326-336)
    # -------------------------------------------------------------------------
    Ift_sq = Vector{Any}(undef, nbranch)
    Itf_sq = Vector{Any}(undef, nbranch)

    for k in 1:nbranch
        f  = f_bus_idx[k]
        t  = t_bus_idx[k]
        θd = theta[f] - theta[t]
        τ  = tau_br[k]

        # Squared current from->to (Dan's MATLAB lines 326-329)
        Ift_sq[k] = @expression(model,
            V[f]^2/τ^4 * (gl[k]^2 + (bl[k]+bsl[k]/2)^2) +
            V[t]^2/τ^2 * (gl[k]^2 + bl[k]^2) -
            2*V[f]*V[t]/τ^3 * (
                (gl[k]^2 + bl[k]*(bl[k]+bsl[k]/2)) * cos(θd - theta_sh[k]) -
                gl[k]*bsl[k]/2                      * sin(θd - theta_sh[k])))

        # Squared current to->from (Dan's MATLAB lines 330-332)
        Itf_sq[k] = @expression(model,
            V[t]^2     * (gl[k]^2 + (bl[k]+bsl[k]/2)^2) +
            V[f]^2/τ^2 * (gl[k]^2 + bl[k]^2) -
            2*V[f]*V[t]/τ * (
                (gl[k]^2 + bl[k]*(bl[k]+bsl[k]/2)) * cos(θd - theta_sh[k]) +
                gl[k]*bsl[k]/2                      * sin(θd - theta_sh[k])))
    end

    # Dan's MATLAB: minimizes obj = (-wft)'*Ift_sq + (-wtf)'*Itf_sq + (-wv)'*V
    # Equivalent here: maximize wft'*Ift_sq + wtf'*Itf_sq + wv'*V
    # res.f = double(-obj) = maximized objective value
    @objective(model, Max,
        sum(wft[k] * Ift_sq[k] for k in 1:nbranch) +
        sum(wtf[k] * Itf_sq[k] for k in 1:nbranch) +
        sum(wv[i]  * V[i]      for i in 1:nbus))

    # -------------------------------------------------------------------------
    # Solve
    # -------------------------------------------------------------------------
    optimize!(model)

    status  = termination_status(model)
    success = status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL

    # -------------------------------------------------------------------------
    # Prepare outputs (matching Dan's MATLAB res structure, lines 368-387)
    # -------------------------------------------------------------------------
    res = Dict{String,Any}(
        "success"            => success,
        "termination_status" => status,
        "f"   => success ? objective_value(model)           : NaN,
        "V"   => success ? value.(V)                        : fill(NaN, nbus),
        "theta" => success ? value.(theta)                  : fill(NaN, nbus),
        "Pg"  => success ? value.(Pg)                       : fill(NaN, ngen),
        "Qg"  => success ? value.(Qg)                       : fill(NaN, ngen),
        "Delta" => success ? value(Delta)                   : NaN,
        "p_ft"  => success ? [value(p_ft[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "q_ft"  => success ? [value(q_ft[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "p_tf"  => success ? [value(p_tf[k]) for k in 1:nbranch] : fill(NaN, nbranch),
        "q_tf"  => success ? [value(q_tf[k]) for k in 1:nbranch] : fill(NaN, nbranch),
    )

    return res
end
