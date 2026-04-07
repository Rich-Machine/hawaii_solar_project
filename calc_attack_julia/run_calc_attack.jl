## Calculate the worst-case line flows and voltage magnitudes for a given set of compromised generators
# Julia port of run_calc_attack.m — direct script translation

using PowerModels
using JuMP
using Ipopt
using Printf
using Plots

include("calc_attack.jl")

verbose = false  # Display solver output?

###############################################################################
## Load test case
###############################################################################
data = PowerModels.parse_file("hawaii40_with_solar.m")

# Remove generators that are turned off
for gk in collect(keys(data["gen"]))
    if get(data["gen"][gk], "gen_status", 1) == 0
        delete!(data["gen"], gk)
    end
end

# Get ordered generator IDs (matching MATPOWER internal ordering after ext2int)
gen_ids    = sort([parse(Int, k) for k in keys(data["gen"])])
branch_ids = sort([parse(Int, k) for k in keys(data["branch"])])
bus_ids    = sort([parse(Int, k) for k in keys(data["bus"])])

nbus    = length(data["bus"])
ngen    = length(data["gen"])
nbranch = length(data["branch"])

###############################################################################
## Choose compromised generators
# MATLAB: c(end-3:end) = true  ->  last 4 generators are compromised
###############################################################################
c = falses(ngen)
c[end-3:end] .= true

# Set the lower generation bounds of the compromised generators to zero
# (models the ability to shut off the generators)
for (idx, id) in enumerate(gen_ids)
    if c[idx]
        data["gen"][string(id)]["pmin"] = 0.0
    end
end

###############################################################################
## Make test case Qd reasonable
# Set reactive power demand based on assumed power factor (matching MATLAB pf=0.96)
###############################################################################
pf = 0.96
for (_, load) in data["load"]
    pd = get(load, "pd", 0.0)
    load["qd"] = pd * sqrt((1/pf^2) - 1)
end

###############################################################################
## Run OPF to get nominal values (prior to attack)
###############################################################################
println("Running baseline OPF...")
res0 = PowerModels.solve_opf(data, ACPPowerModel, Ipopt.Optimizer)

# Update data with baseline OPF solution
for (bk, bsol) in res0["solution"]["bus"]
    haskey(data["bus"], bk) || continue
    data["bus"][bk]["vm"] = bsol["vm"]
    data["bus"][bk]["va"] = bsol["va"]
end
for (gk, gsol) in res0["solution"]["gen"]
    haskey(data["gen"], gk) || continue
    data["gen"][gk]["pg"] = gsol["pg"]
    data["gen"][gk]["qg"] = gsol["qg"]
    gen_bus = data["gen"][gk]["gen_bus"]
    bk = string(gen_bus)
    if haskey(res0["solution"]["bus"], bk)
        data["gen"][gk]["vg"] = res0["solution"]["bus"][bk]["vm"]
    end
end

###############################################################################
## Set S-curve parameters and participation factors
###############################################################################
pvpq_k = 100.0  # Slope parameter for smoothed PV/PQ switching characteristic

# Participation factors based on nominal active power outputs
Pg0    = [data["gen"][string(id)]["pg"] for id in gen_ids]
Pg_tot = sum(Pg0)
alpha  = Pg_tot > 0 ? Pg0 ./ Pg_tot : zeros(ngen)

###############################################################################
## Find worst-case outcomes for the given set of compromised generators
###############################################################################

# Branch thermal ratings for percentage-loading calculation
rate_a = [get(data["branch"][string(id)], "rate_a", Inf) for id in branch_ids]

# Placeholder variables to store results
Ift  = fill(NaN, nbranch)
Itf  = fill(NaN, nbranch)
Vlow = fill(NaN, nbus)
Vhigh= fill(NaN, nbus)

## Maximize flows in f->t direction (commented out, matching MATLAB)
for i in 1:nbranch
    wft = zeros(nbranch); wtf = zeros(nbranch); wv = zeros(nbus)
    wft[i] = 1.0
    @printf("Computing Worst-Case Flow (from -> to direction) for line %d of %d\n", i, nbranch)
    res = calc_attack(data, wft, wtf, wv, c, alpha, pvpq_k; verbose=verbose)
    if res["success"]
        Ift[i] = sqrt(res["f"]) * 100 / rate_a[i]
    end
end

## Maximize flows in t->f direction (commented out, matching MATLAB)
for i in 1:nbranch
    wft = zeros(nbranch); wtf = zeros(nbranch); wv = zeros(nbus)
    wtf[i] = 1.0
    @printf("Computing Worst-Case Flow (to -> from direction) for line %d of %d\n", i, nbranch)
    res = calc_attack(data, wft, wtf, wv, c, alpha, pvpq_k; verbose=verbose)
    if res["success"]
        Itf[i] = sqrt(res["f"]) * 100 / rate_a[i]
    end
end

max_violations = max.(Ift, Itf)

## Minimize voltages
for i in 1:nbus
    wft = zeros(nbranch); wtf = zeros(nbranch); wv = zeros(nbus)
    wv[i] = -1.0
    @printf("Computing Worst-Case Lower Voltage for bus %d of %d\n", i, nbus)
    res = calc_attack(data, wft, wtf, wv, c, alpha, pvpq_k; verbose=verbose)
    if res["success"]
        Vlow[i] = -res["f"]
        println(Vlow)
    end
end

## Maximize voltages (commented out, matching MATLAB)
for i in 1:nbus
    wft = zeros(nbranch); wtf = zeros(nbranch); wv = zeros(nbus)
    wv[i] = 1.0
    @printf("Computing Worst-Case Upper Voltage for bus %d of %d\n", i, nbus)
    res = calc_attack(data, wft, wtf, wv, c, alpha, pvpq_k; verbose=verbose)
    if res["success"]
        Vhigh[i] = res["f"]
    end
end

###############################################################################
## Plotting — helper for shaded vertical bands (equivalent to MATLAB patch)
###############################################################################

# Add a shaded rectangle between x1 and x2 spanning full plot height.
# Call this AFTER the histogram so ylims are set.
function vband!(p, x1, x2, color; fillalpha=0.12, label="")
    yl_lo, yl_hi = ylims(p)
    plot!(p, Shape([x1, x2, x2, x1], [yl_lo, yl_lo, yl_hi, yl_hi]);
          fillcolor=color, fillalpha=fillalpha,
          linewidth=0, linecolor=:transparent, label=label)
    ylims!(p, yl_lo, yl_hi)  # Restore y limits so the shape doesn't expand the axes
end

###############################################################################
## Plotting — Line violations
# Requires the line flow loops above to be uncommented.
###############################################################################
if !all(isnan.(Ift)) || !all(isnan.(Itf))
    max_violations = max.(Ift, Itf)
    valid_mv = filter(!isnan, max_violations)
    xl_lo = minimum(valid_mv)
    xl_hi = maximum(valid_mv)

    pL = histogram(max_violations; bins=20,
                   color=:gray, fillalpha=0.7, label="",
                   title="Worst-Case Line Loading", titlefontsize=14,
                   xlabel="Percentage loading (%)", ylabel="Number of branches",
                   labelfontsize=14, grid=true, legend=:topright)

    vband!(pL, xl_lo, 100.0,  :blue; label="Normal operation")
    vband!(pL, 100.0,  xl_hi,  :red;  label="Constraint violation")
    vline!([100.0]; color=:red, linestyle=:dash, linewidth=2, label="")

    display(pL)
end

###############################################################################
## Plotting — Voltage violations (two stacked subplots)
###############################################################################
valid_vhigh = filter(!isnan, Vhigh)
valid_vlow  = filter(!isnan, Vlow)

# Top panel: Overvoltage (Vhigh)
# Uncomment the Vhigh loop above to populate this panel.
if isempty(valid_vhigh)
    p1 = plot(; title="Worst-Case Overvoltage\n(uncomment Vhigh loop to run)",
               titlefontsize=14, grid=true)
else
    x1_hi = maximum(valid_vhigh)
    p1 = histogram(Vhigh; bins=20,
                   color=:gray, fillalpha=0.7, label="",
                   title="Worst-Case Overvoltage", titlefontsize=16,
                   xlabel="Voltage (p.u.)", ylabel="Number of buses",
                   labelfontsize=14, grid=true, legend=:topleft)

    vband!(p1, 1.00, 1.05,   RGB(0.3, 0.6, 1.0); label="Normal operation")
    vband!(p1, 1.05, x1_hi,  RGB(1.0, 0.4, 0.4); label="Constraint violation")
    vline!([1.05]; color=:red,  linestyle=:dash, linewidth=2, label="")
    vline!([1.00]; color=:blue, linewidth=1.2,               label="")
end

# Bottom panel: Undervoltage (Vlow)
if isempty(valid_vlow)
    p2 = plot(; title="Worst-Case Undervoltage\n(uncomment Vlow loop to run)",
               titlefontsize=14, grid=true)
else
    x2_lo = minimum(valid_vlow)
    p2 = histogram(Vlow; bins=20,
                   color=:gray, fillalpha=0.7, label="",
                   title="Worst-Case Undervoltage", titlefontsize=16,
                   xlabel="Voltage (p.u.)", ylabel="Number of buses",
                   labelfontsize=14, grid=true, legend=:topright)

    vband!(p2, x2_lo, 0.95,  RGB(1.0, 0.4, 0.4); label="Constraint violation")
    vband!(p2, 0.95,  1.00,  RGB(0.3, 0.6, 1.0); label="Normal operation")
    vline!([0.95]; color=:red,  linestyle=:dash, linewidth=2, label="")
    vline!([1.00]; color=:blue, linewidth=1.2,               label="")
end

p_volt = plot(p1, p2;
              layout=(2, 1), size=(800, 800),
              plot_title="Worst-Case Voltage Analysis — Distributions",
              plot_titlefontsize=18)
display(p_volt)
