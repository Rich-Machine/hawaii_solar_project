using PowerModels
using JuMP
using Ipopt
using CSV
using DataFrames
using PowerPlots

# network_data = PowerModels.parse_file("hawaii40.m")

# data = CSV.File("oahu_solar_inverters.csv") |> DataFrame


function add_solar_inverters!(network_data, inverter_data)
    active_inverters = []
    matched_inverters = []
    for row in eachrow(inverter_data)
        if row[:Capacity] > 0 
            push!(active_inverters, row[:Column1])
            for b in 1:length(network_data["bus"])
                if abs(row[:Latitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[1])) <= 1e-1 && 
                    abs(row[:Longitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[2])) <= 1e-1 && network_data["bus"][string(b)]["base_kv"] == 69.0 && network_data["bus"][string(b)]["bus_type"] == 2 && !(row[:Column1] in matched_inverters)
                    push!(matched_inverters, row[:Column1])
                    println("Matched inverter at (", row[:Column1], ") to bus ", network_data["bus"][string(b)]["name"], "in first condition")
                    gen_id = 100 + row[:Column1]
                    network_data["gen"][string(gen_id)] = Dict(
                        "gen_bus" => b,
                        "pg" => row[:Capacity] * 1e-5,
                        "qg" => row[:Capacity] * 1e-5 * 0.2,
                        "pmax" => row[:Capacity] * 1e-5,
                        "pmin" => 0.0,
                        "qmax" => row[:Capacity] * 1e-5,
                        "qmin" => -row[:Capacity] * 1e-5,
                        "gen_status" => 1,
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

                elseif abs(row[:Latitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[1])) > 1e-1 && abs(row[:Latitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[1])) < 4e-1 &&
                       abs(row[:Longitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[2])) > 1e-1 && abs(row[:Longitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[2])) < 4e-1 && network_data["bus"][string(b)]["base_kv"] == 69.0 && network_data["bus"][string(b)]["bus_type"] == 2 && !(row[:Column1] in matched_inverters)
                    push!(matched_inverters, row[:Column1])
                    println("Matched inverter at (", row[:Column1], ") to bus ", network_data["bus"][string(b)]["name"], " in second condition")
                    gen_id = 100 + row[:Column1]
                    network_data["gen"][string(gen_id)] = Dict(
                        "gen_bus" => b,
                        "pg" => row[:Capacity] * 1e-5,
                        "qg" => row[:Capacity] * 1e-5 * 0.2,
                        "pmax" => row[:Capacity] * 1e-5,
                        "pmin" => 0.0,
                        "qmax" => row[:Capacity] * 1e-5,
                        "qmin" => -row[:Capacity] * 1e-5,
                        "gen_status" => 1,
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

                elseif abs(row[:Latitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[1])) > 4e-1 && abs(row[:Latitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[1])) < 8e-1 &&
                    abs(row[:Longitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[2])) > 4e-1 && abs(row[:Longitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[2])) < 8e-1 && network_data["bus"][string(b)]["base_kv"] == 69.0 && network_data["bus"][string(b)]["bus_type"] == 2 && !(row[:Column1] in matched_inverters)
                    push!(matched_inverters, row[:Column1])
                    println("Matched inverter at (", row[:Column1], ") to bus ", network_data["bus"][string(b)]["name"], " in third condition")
                    gen_id = 100 + row[:Column1]
                    network_data["gen"][string(gen_id)] = Dict(
                        "gen_bus" => b,
                        "pg" => row[:Capacity] * 1e-5,
                        "qg" => row[:Capacity] * 1e-5 * 0.2,
                        "pmax" => row[:Capacity] * 1e-5,
                        "pmin" => 0.0,
                        "qmax" => row[:Capacity] * 1e-5,
                        "qmin" => -row[:Capacity] * 1e-5,
                        "gen_status" => 1,
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



                elseif abs(row[:Latitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[1])) > 8e-1 && 
                    abs(row[:Longitude] - parse(Float64, split(network_data["coordinates"][string(b)]["col_1"], ",")[2])) > 8e-1 && network_data["bus"][string(b)]["base_kv"] == 69.0 && network_data["bus"][string(b)]["bus_type"] == 2 && !(row[:Column1] in matched_inverters)
                    push!(matched_inverters, row[:Column1])
                    println("Matched inverter at (", row[:Column1], ") to bus ", network_data["bus"][string(b)]["name"], " in fourth condition")
                    gen_id = 100 + row[:Column1]
                    network_data["gen"][string(gen_id)] = Dict(
                        "gen_bus" => b,
                        "pg" => row[:Capacity] * 1e-5,
                        "qg" => row[:Capacity] * 1e-5 * 0.2,
                        "pmax" => row[:Capacity] * 1e-5,
                        "pmin" => 0.0,
                        "qmax" => row[:Capacity] * 1e-5,
                        "qmin" => -row[:Capacity] * 1e-5,
                        "gen_status" => 1,
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
    print("Added $(length(matched_inverters)) solar inverters to the network.\n")
    print(sum(row[:Capacity] for row in eachrow(inverter_data) if row[:Column1] in matched_inverters) , " MW of solar capacity added.\n")
end

# add_solar_inverters!(network_data, data)

