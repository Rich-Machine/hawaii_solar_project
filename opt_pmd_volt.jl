using PowerModels
using Ipopt
using JuMP
using PowerPlots
using CSV
using DataFrames
using Gurobi
using Plots
# using LaTeXStrings
const optimizer = Gurobi.Optimizer
include("add_solar_inverters.jl")

## Load the model and solve the base case OPF
eng = PowerModels.parse_file("hawaii40.m")
add_solar_inverters!(eng, CSV.File("oahu_solar_inverters.csv") |> DataFrame)
results = solve_opf(eng, ACPPowerModel, Ipopt.Optimizer)

## Plot the network
powerplot(eng, height=500, width=500, bus_size=5, gen_size=75, load_size=0.1, branch_size=2)

for b in keys(eng["bus"])
    eng["bus"][b]["vmax"] = 1.2
    eng["bus"][b]["vmin"] = 0.8
end

solar = []
fuel = []
for s in keys(eng["gen"])
    if eng["genfuel"][s]["col_1"] == "solar"
        push!(solar, parse(Int, s))
    elseif s == "41" || s == "33" || s == "43" || s == "40" || s == "30" || s == "45"
        continue
    else
        push!(fuel, parse(Int, s))
        println("Fuel Generator: $s")
    end
end

gen = collect(1:length(eng["gen"]))
buses = collect(1:length(eng["bus"]))
no_gens = length(eng["gen"])

budget = 0.2
delta = 0.01

deltas = [0.01]
budgets = [0.5]
no_of_violations = []

for budget in budgets
    for delta in deltas
        pm = instantiate_model(eng, ACPPowerModel, PowerModels.build_opf)
        opt_model = pm.model

        ## Add the solar panel attack constraints; three phases per solar
        begin
            JuMP.@variable(opt_model, u_s[s in solar], binary=true)
            JuMP.@variable(opt_model, P_s[s in solar])
            JuMP.@variable(opt_model, Q_s[s in solar])
            JuMP.@constraint(opt_model, [s in solar], (1 - u_s[s]) * eng["gen"]["$s"]["pmax"][1] <= P_s[s])
            JuMP.@constraint(opt_model, [s in solar], P_s[s] <= eng["gen"]["$s"]["pmax"][1])
            JuMP.@constraint(opt_model, [s in solar], eng["gen"]["$s"]["qmin"][1] * u_s[s] <= Q_s[s])
            JuMP.@constraint(opt_model, [s in solar], Q_s[s] <= eng["gen"]["$s"]["qmax"][1] * u_s[s])
            # JuMP.@constraint(opt_model, [s in solar], (1 - u_s[s]) * results["solution"]["gen"]["$s"]["pg"] <= P_s[s])
            # JuMP.@constraint(opt_model, [s in solar], P_s[s] <= eng["gen"]["$s"]["pmax"][1] * u_s[s] + results["solution"]["gen"]["$s"]["pg"] * (1 - u_s[s]))
            # JuMP.@constraint(opt_model, [s in solar], results["solution"]["gen"]["$s"]["qg"] * (1 - u_s[s]) + eng["gen"]["$s"]["qmin"][1] * u_s[s] <= Q_s[s])
            # JuMP.@constraint(opt_model, [s in solar], Q_s[s] <= eng["gen"]["$s"]["qmax"][1] * u_s[s] + results["solution"]["gen"]["$s"]["qg"] * (1 - u_s[s]))

        end

        begin
            JuMP.@variable(opt_model, voltage[b in buses])
            JuMP.@variable(opt_model, y_over[b in buses], Bin)
            JuMP.@variable(opt_model, y_under[b in buses], Bin)
            JuMP.@variable(opt_model, P_g[g in fuel])
            JuMP.@variable(opt_model, Q_g[g in fuel])
            JuMP.@constraint(opt_model, [g in fuel], P_g[g] == results["solution"]["gen"]["$g"]["pg"])
            JuMP.@constraint(opt_model, [g in fuel], Q_g[g] == results["solution"]["gen"]["$g"]["qg"])

            JuMP.@constraint(opt_model, [b in buses], voltage[b] == var(pm, :vm)[b])
            JuMP.@constraint(opt_model, [b in buses], voltage[b] >= (1.1 + delta) * y_over[b])
            JuMP.@constraint(opt_model, [b in buses], voltage[b] <= (0.9 - delta) * y_under[b] + 1000(1 - y_under[b]))
            JuMP.@constraint(opt_model, [b in buses],  y_over[b] + y_under[b] <= 1)
            JuMP.@constraint(opt_model, [s in solar], sum(u_s) <= budget * length(solar))  
            JuMP.@constraint(opt_model, [s in solar], var(pm, :pg, s) .== P_s[s])   
            JuMP.@constraint(opt_model, [s in solar], var(pm, :qg, s) .== Q_s[s])   
            JuMP.@constraint(opt_model, [g in fuel], var(pm, :pg, g) .== P_g[g])   
            JuMP.@constraint(opt_model, [g in fuel], var(pm, :qg, g) .== Q_g[g])   

            JuMP.@objective(opt_model, Max, sum((y_over[b]) for b in buses) + sum((y_under[b]) for b in buses))
        end

        ## Solve the model
        begin
            result = optimize_model!(pm, optimizer=Gurobi.Optimizer)
            # println(JuMP.value.(u_s))
            # println(JuMP.value.(y_over))
            # println(JuMP.value.(y_under))
            print
            println(round.(JuMP.value.(voltage), digits=4))
            println("Objective value: ", JuMP.objective_value(opt_model))
            println("Total violations: ", sum(JuMP.value.(y_over)) + sum(JuMP.value.(y_under)), " out of ", length(buses))
            println("Total attacked solar: ", sum(JuMP.value.(u_s)), " out of ", length(solar))
            println("\n")
        end

        for s in solar
            # println("Solar Panel $s attacked: ", JuMP.value(u_s[s]), " Pmax: ", round(eng["gen"]["$s"]["pmax"][1], digits=4) )
            println("Solar Panel $s attacked: ", JuMP.value(u_s[s]), " Pg: ", round(JuMP.value(P_s[s]), digits=4), " Qg: ", round(JuMP.value(Q_s[s]), digits=4) )
        end
        println("\n")
        push!(no_of_violations, JuMP.objective_value(opt_model))
    end
end

no_of_violations = reshape(no_of_violations, (length(deltas), length(budgets)))

# Plotting the results
# plot(budgets, no_of_violations[:, 2], label=L"\Delta = 0.02", xlabel="Budget", ylabel="Number of Violations", title="Voltage Violations vs Budget", legend=:bottomright, marker=:circle, linewidth=2)
# plot!(budgets, no_of_violations[:, 4], label=L"\Delta = 0.04", marker=:asterisk, linewidth=2)
# plot!(budgets, no_of_violations[:, 6], label=L"\Delta = 0.06", marker=:square, linewidth=2)
# plot!(budgets, no_of_violations[:, 8], label=L"\Delta = 0.08", marker=:diamond, linewidth=2)
# plot!(budgets, no_of_violations[:, 10], label=L"\Delta = 0.10", marker=:star5, linewidth=2)

# plot(budgets, no_of_violations[1, :], label=L"\Delta = 0.01", xlabel="Budget", ylabel="Number of Violations", title="Voltage Violations vs Budget", legend=:bottomright, marker=:circle, linewidth=2)
# plot(deltas, no_of_violations[:, 1], label="Budget = 20%", xlabel=L"\textbf{Delta,} \Delta", ylabel="Number of Violations", title="Voltage Violations vs Delta", legend=:topright, marker=:circle, linewidth=2)
# plot!(deltas, no_of_violations[:, 2], label="Budget = 40%", marker=:asterisk, linewidth=2)
# plot!(deltas, no_of_violations[:, 3], label="Budget = 60%", marker=:square, linewidth=2)
# plot!(deltas, no_of_violations[:, 4], label="Budget = 80%", marker=:diamond, linewidth=2)
# plot!(deltas, no_of_violations[:, 5], label="Budget = 100%", marker=:star5, linewidth=2)

# savefig("voltage_violations_vs_budget_few.png")

for g in fuel
    if "g" in keys(results["solution"]["gen"])
        continue
    else
        println("Generator $g not in the solution.")
    end
end

