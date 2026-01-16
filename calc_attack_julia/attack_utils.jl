"""
    attack_utils.jl

Utility functions for the calc_attack Julia port:
- S-curve parameter fitting for AGC model
"""

"""
    optimal_scurve_k_asym_brute_force(Pgmin, Pgmax, Pg0, alpha)

Compute optimal S-curve parameters (k, Δ0) that minimize L1 error with respect
to the piecewise linear AGC characteristic.

# Arguments
- `Pgmin`: Minimum active power output
- `Pgmax`: Maximum active power output
- `Pg0`: Nominal active power output
- `alpha`: Participation factor (slope of linear AGC region)

# Returns
- `k_star`: Optimal S-curve slope parameter
- `Delta0_star`: Optimal S-curve center offset

The S-curve model is: Pg = Pgmin + (Pgmax - Pgmin) / (1 + exp(-k*(Δ - Δ0)))
"""
function optimal_scurve_k_asym_brute_force(Pgmin::Float64, Pgmax::Float64,
                                           Pg0::Float64, alpha::Float64)
    # Grid of Delta values for error evaluation
    Delta = collect(-2.0:0.01:2.0)

    # Grid of k values, scaled by alpha/(Pgmax - Pgmin)
    k_base = collect(0.001:0.001:10.0)
    k = k_base .* alpha / (Pgmax - Pgmin)

    err = fill(NaN, length(k))

    # Piecewise linear characteristic breakpoints
    Delta1 = (Pgmin - Pg0) / alpha  # Where saturation begins (low)
    Delta2 = (Pgmax - Pg0) / alpha  # Where saturation begins (high)

    for i in 1:length(k)
        # Compute S-curve center for this k value
        # This ensures the S-curve passes through (0, Pg0)
        if Pg0 <= Pgmin || Pg0 >= Pgmax
            # Degenerate case: nominal output at limits
            err[i] = Inf
            continue
        end
        Delta0 = (1/k[i]) * log((Pgmax - Pg0) / (Pg0 - Pgmin))

        # Compute S-curve values
        Pg = Pgmin .+ (Pgmax - Pgmin) ./ (1.0 .+ exp.(-k[i] .* (Delta .- Delta0)))

        # Compute error vs piecewise linear characteristic
        # Region 1: Delta < Delta1 -> Pg should be Pgmin
        # Region 2: Delta1 <= Delta <= Delta2 -> Pg should be Pg0 + alpha*Delta
        # Region 3: Delta > Delta2 -> Pg should be Pgmax

        total_err = 0.0
        for j in 1:length(Delta)
            if Delta[j] < Delta1
                total_err += abs(Pg[j] - Pgmin)
            elseif Delta[j] > Delta2
                total_err += abs(Pg[j] - Pgmax)
            else
                total_err += abs(Pg[j] - (Pg0 + alpha * Delta[j]))
            end
        end
        err[i] = total_err
    end

    # Find minimum error
    _, idx = findmin(err)
    k_star = k[idx]
    Delta0_star = (1/k_star) * log((Pgmax - Pg0) / (Pg0 - Pgmin))

    return k_star, Delta0_star
end
