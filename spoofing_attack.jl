using PowerModels
using Ipopt
using JuMP
using PowerPlots
using CSV
using DataFrames
using Gurobi
using Plots
using StatsPlots
using Random
# Random.seed!(1234)

# using LaTeXStrings
const optimizer = Gurobi.Optimizer
include("add_solar_inverters.jl")

## Load the model and solve the base case OPF
eng = PowerModels.parse_file("hawaii40.m")
add_solar_inverters!(eng, CSV.File("oahu_solar_inverters.csv") |> DataFrame)
base_results = solve_opf(eng, ACPPowerModel, Ipopt.Optimizer)

# solar = []
# for s in keys(eng["gen"])
#     if eng["genfuel"][s]["col_1"] == "solar" && parse(Int, s) > 100
#         push!(solar, parse(Int, s))
#     end
# end

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

sum(eng["gen"]["$s"]["pmax"] for s in solar)

# # scales = [0, 1, 10, 100, 1000]
# scales = [0, 1]
# vm_matrix = []
# va_matrix = []
# for scale in scales
#     if scale == 0
#         continue
#         for s in solar
#             display(eng["gen"]["$s"]["pg"])
#             display(eng["gen"]["$s"]["qg"])
#         end
#     else
#         for s in solar
#             eng["gen"]["$s"]["pmax"] = eng["gen"]["$s"]["pmax"] .* [scale]
#             eng["gen"]["$s"]["pg"] = rand() * eng["gen"]["$s"]["pmax"]
#             eng["gen"]["$s"]["qg"] = rand() * eng["gen"]["$s"]["qmax"]
#             # eng["gen"]["$s"]["qg"] = rand() * (eng["gen"]["$s"]["qmax"][1] - eng["gen"]["$s"]["qmin"][1]) + eng["gen"]["$s"]["qmin"][1]
#         end
#     end

#     results = solve_pf(eng, ACPPowerModel, Ipopt.Optimizer)
#     voltage_magnitudes = []
#     voltage_angles = []
#     push!(voltage_magnitudes, [results["solution"]["bus"][string(b)]["vm"] for b in keys(eng["bus"])]...)
#     push!(voltage_angles, [results["solution"]["bus"][string(b)]["va"] for b in keys(eng["bus"])]...)
#     voltage_magnitudes = reshape(voltage_magnitudes, (length(eng["bus"]), 1))
#     voltage_angles = reshape(voltage_angles, (length(eng["bus"]), 1))
#     append!(vm_matrix, voltage_magnitudes)
#     append!(va_matrix, voltage_angles)
# end

# vm_matrix = reshape(vm_matrix, (length(eng["bus"]), length(scales)))
# va_matrix = reshape(va_matrix, (length(eng["bus"]), length(scales)))

# StatsPlots.boxplot(vm_matrix[:,1], ylabel="Voltage magnitude in p.u.", outliers=false)
# StatsPlots.boxplot!(vm_matrix[:,2],  outliers=false)
# # StatsPlots.boxplot(va_matrix[:,1] * 57.2957795, xlabel="Under attack", outliers=false)


# scales = [0, 1, 10, 100, 1000]
scales = [0]
vm_matrix = []
va_matrix = []
for scale in scales
    # for s in solar
    #     eng["gen"]["$s"]["pmax"] = eng["gen"]["$s"]["pmax"] .* [0]
    #     eng["gen"]["$s"]["pg"] = base_results["solution"]["gen"]["$s"]["pg"]
    #     eng["gen"]["$s"]["qg"] = base_results["solution"]["gen"]["$s"]["pg"]
    # end



    results = solve_pf(eng, ACPPowerModel, Ipopt.Optimizer)
    voltage_magnitudes = []
    voltage_angles = []
    push!(voltage_magnitudes, [results["solution"]["bus"][string(b)]["vm"] for b in keys(eng["bus"])]...)
    push!(voltage_angles, [results["solution"]["bus"][string(b)]["va"] for b in keys(eng["bus"])]...)
    voltage_magnitudes = reshape(voltage_magnitudes, (length(eng["bus"]), 1))
    voltage_angles = reshape(voltage_angles, (length(eng["bus"]), 1))
    append!(vm_matrix, voltage_magnitudes)
    append!(va_matrix, voltage_angles)
end

vm_matrix = reshape(vm_matrix, (length(eng["bus"]), length(scales)))
va_matrix = reshape(va_matrix, (length(eng["bus"]), length(scales)))

StatsPlots.boxplot(vm_matrix[:,1], ylabel="Voltage magnitude", outliers=false, title=" Nodal Voltage Magnitude under the different scenarios", label="Base Case")
# StatsPlots.boxplot(va_matrix[:,1] * 57.2957795, ylabel="Voltage angle", outliers=false, title=" Nodal Voltage Angles under the different scenarios", label="Base Case")


# scales = [0, 1, 10, 100, 1000]
scales = [1]
vm_matrix = []
va_matrix = []
for scale in scales
    for s in solar
        eng["gen"]["$s"]["pmax"] = eng["gen"]["$s"]["pmax"] .* [scale]
        eng["gen"]["$s"]["pg"] = rand() * eng["gen"]["$s"]["pmax"][1] * rand((-1, 1))
        eng["gen"]["$s"]["qg"] = rand() * eng["gen"]["$s"]["qmax"][1] * rand((-1, 1))
    end
    # for g in fuel
    #     eng["gen"]["$g"]["pg"] = base_results["solution"]["gen"]["$g"]["pg"]
    #     eng["gen"]["$g"]["qg"] = base_results["solution"]["gen"]["$g"]["qg"]
    # end

    results = solve_pf(eng, ACPPowerModel, Ipopt.Optimizer)
    display(results["solution"]["gen"]["1"]["pg"])

    voltage_magnitudes = []
    voltage_angles = []
    push!(voltage_magnitudes, [results["solution"]["bus"][string(b)]["vm"] for b in keys(eng["bus"])]...)
    push!(voltage_angles, [results["solution"]["bus"][string(b)]["va"] for b in keys(eng["bus"])]...)
    voltage_magnitudes = reshape(voltage_magnitudes, (length(eng["bus"]), 1))
    voltage_angles = reshape(voltage_angles, (length(eng["bus"]), 1))
    append!(vm_matrix, voltage_magnitudes)
    append!(va_matrix, voltage_angles)
end

vm_matrix = reshape(vm_matrix, (length(eng["bus"]), length(scales)))
va_matrix = reshape(va_matrix, (length(eng["bus"]), length(scales)))

StatsPlots.boxplot!(vm_matrix[:,1], outliers=false, label="1x Solar Penetration")
# StatsPlots.boxplot!(va_matrix[:,1] * 57.2957795, ylabel="Voltage angles", outliers=false, label="1x Solar Penetration")


## Load the model and solve the base case OPF
eng = PowerModels.parse_file("hawaii40.m")
add_solar_inverters!(eng, CSV.File("oahu_solar_inverters.csv") |> DataFrame)

solar = []
for s in keys(eng["gen"])
    if eng["genfuel"][s]["col_1"] == "solar"
        push!(solar, parse(Int, s))
    end
end

# scales = [0, 1, 10, 100, 1000]
scales = [10]
vm_matrix = []
va_matrix = []
for scale in scales
    for s in solar
        eng["gen"]["$s"]["pmax"] = eng["gen"]["$s"]["pmax"] .* [scale]
        eng["gen"]["$s"]["pg"] = rand() * eng["gen"]["$s"]["pmax"][1]* rand((-1, 1))
        eng["gen"]["$s"]["qg"] = rand() * eng["gen"]["$s"]["qmax"][1]* rand((-1, 1))
    end

    # for g in fuel
    #     eng["gen"]["$g"]["pg"] = base_results["solution"]["gen"]["$g"]["pg"]
    #     eng["gen"]["$g"]["qg"] = base_results["solution"]["gen"]["$g"]["qg"]
    # end

    results = solve_pf(eng, ACPPowerModel, Ipopt.Optimizer)
    display(results["solution"]["gen"]["1"]["pg"])

    voltage_magnitudes = []
    voltage_angles = []
    push!(voltage_magnitudes, [results["solution"]["bus"][string(b)]["vm"] for b in keys(eng["bus"])]...)
    push!(voltage_angles, [results["solution"]["bus"][string(b)]["va"] for b in keys(eng["bus"])]...)
    voltage_magnitudes = reshape(voltage_magnitudes, (length(eng["bus"]), 1))
    voltage_angles = reshape(voltage_angles, (length(eng["bus"]), 1))
    append!(vm_matrix, voltage_magnitudes)
    append!(va_matrix, voltage_angles)
end

vm_matrix = reshape(vm_matrix, (length(eng["bus"]), length(scales)))
va_matrix = reshape(va_matrix, (length(eng["bus"]), length(scales)))

StatsPlots.boxplot!(vm_matrix[:,1], outliers=false, label="10x Solar Penetration")
# StatsPlots.boxplot!(va_matrix[:,1] * 57.2957795, outliers=false, label="10x Solar Penetration")





## Load the model and solve the base case OPF
eng = PowerModels.parse_file("hawaii40.m")
add_solar_inverters!(eng, CSV.File("oahu_solar_inverters.csv") |> DataFrame)
# base_results = solve_opf(eng, ACPPowerModel, Ipopt.Optimizer)

solar = []
for s in keys(eng["gen"])
    if eng["genfuel"][s]["col_1"] == "solar"
        push!(solar, parse(Int, s))
    end
end

# scales = [0, 1, 10, 100, 1000]
scales = [20]
vm_matrix = []
va_matrix = []
for scale in scales
    for s in solar
        eng["gen"]["$s"]["pmax"] = eng["gen"]["$s"]["pmax"] .* [scale]
        eng["gen"]["$s"]["pg"] = rand() * eng["gen"]["$s"]["pmax"][1] * rand((-1, 1))
        eng["gen"]["$s"]["qg"] = rand() * eng["gen"]["$s"]["qmax"][1] * rand((-1, 1))
    end

    # for g in fuel
    #     eng["gen"]["$g"]["pg"] = base_results["solution"]["gen"]["$g"]["pg"]
    #     eng["gen"]["$g"]["qg"] = base_results["solution"]["gen"]["$g"]["qg"]
    # end

    results = solve_pf(eng, ACPPowerModel, Ipopt.Optimizer)
    display(results["solution"]["gen"]["1"]["pg"])
    voltage_magnitudes = []
    voltage_angles = []
    push!(voltage_magnitudes, [results["solution"]["bus"][string(b)]["vm"] for b in keys(eng["bus"])]...)
    push!(voltage_angles, [results["solution"]["bus"][string(b)]["va"] for b in keys(eng["bus"])]...)
    voltage_magnitudes = reshape(voltage_magnitudes, (length(eng["bus"]), 1))
    voltage_angles = reshape(voltage_angles, (length(eng["bus"]), 1))
    append!(vm_matrix, voltage_magnitudes)
    append!(va_matrix, voltage_angles)
end

vm_matrix = reshape(vm_matrix, (length(eng["bus"]), length(scales)))
va_matrix = reshape(va_matrix, (length(eng["bus"]), length(scales)))
StatsPlots.boxplot!(vm_matrix[:,1], outliers=false, label="100x Solar Penetration")
hline!([0.95], linestyle=:dash, color=:red, label="Voltage Magnitude Limit")
# StatsPlots.boxplot!(va_matrix[:,1] * 57.2957795, outliers=false, label="20x Solar Penetration")
# hline!([30], linestyle=:dash, color=:red, label="Voltage Angle Limit")



# savefig("voltage_angle_spoofing_attack.pdf")


# eng = PowerModels.parse_file("hawaii40.m")
# add_solar_inverters!(eng, CSV.File("oahu_solar_inverters.csv") |> DataFrame)
# sum(eng["gen"]["$s"]["pmax"] for s in solar)
# sum(eng["gen"]["$s"]["pg"] for s in fuel)
# sum(eng["load"]["$b"]["pd"] for b in keys(eng["load"]))

# base_results = solve_opf(eng, ACPPowerModel, Ipopt.Optimizer)