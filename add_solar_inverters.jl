using PowerModels
using JuMP
using Ipopt
using CSV
using DataFrames
using PowerPlots

network_data = PowerModels.parse_file("hawaii40.m")

# solve_opf(network_data, ACPPowerModel, Ipopt.Optimizer)
solve_opf(network_data, DCPPowerModel, Ipopt.Optimizer)



data = CSV.File("oahu_solar_inverters.csv") |> DataFrame

function add_solar_inverters!(network_data, inverter_data)
    for row in eachrow(inverter_data)
        if row[:Capacity] > 10
            for b in 1:length(network_data["bus"])
                if abs(row[:Latitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[1])) < 1e-1 && 
                    abs(row[:Longitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[2])) < 1e-1 && network_data["bus"][string(b)]["base_kv"] == 69.0 && network_data["bus"][string(b)]["bus_type"] == 2
                    # row[:bus_id] = b
                    println("Matched inverter at (", row[:Latitude], ", ", row[:Longitude], ") to bus ", network_data["bus"][string(b)]["name"])
                    # gen_id = "solar_inverter_" * string(row[:Column1])
                    gen_id = 100 + row[:Column1]
                    network_data["gen"][string(gen_id)] = Dict(
                        "gen_bus" => b,
                        "pg" => row[:Capacity] * 10e-5,
                        "qg" => row[:Capacity] * 10e-5 * 0.2,
                        "pmax" => row[:Capacity] * 10e-5,
                        "pmin" => 0.0,
                        # "qmax" => row[:Capacity] * 10e-5 * 0.2,
                        # "qmin" => -row[:Capacity] * 10e-5 * 0.2,
                        "qmax" => row[:Capacity] * 10e-5,
                        "qmin" => -row[:Capacity] * 10e-5,
                        "gen_status" => 1,
                        # "gen_type" => "solar_inverter",
                        "ncost" => 0,
                        "vg"  => 1.0,
                        "cost" => Float64[],
                        "qc1max"     => 0.0,
                        "model"      => 2,
                        "shutdown"   => 0.0,
                        "startup"    => 0.0,
                        "qc2max"     => 0.0,
                        "ramp_agc"   => 0.0,
                        "ramp_10"    => 0.0,
                        "mbase"      => 100 ,
                        "source_id"  => Any["gen", row[:Column1]],
                        "mu_pmax"    => 0.0,
                        "pc2"        => 0.0,
                        "mu_pmin"    => 0.0,
                        "index"      => 100 + row[:Column1],
                        "qc1min"     => 0.0,
                        "qc2min"     => 0.0,
                        "pc1"        => 0.0,
                        "ramp_q"     => 0.0,
                        "mu_qmax"    => 0.0,
                        "ramp_30"    => 0.0,
                        "mu_qmin"    => 0.0,
                        "apf"        => 0.0
                    )
                    network_data["genfuel"][string(gen_id)] = Dict(
                        "source_id" => Any["genfuel", gen_id],
                        "col_1"     => "solar",
                        "index"     => gen_id
                    )
                    continue

                    return network_data
                end
            end
        end
    end
end


add_solar_inverters!(network_data, data)