# Power-System Analysis for *Grid Trouble in Paradise*

This repository contains the power-system analysis code accompanying the
paper. It takes the Oahu-localized exposed-inverter set produced by the
measurement pipeline, places those inverters on a synthetic Oahu transmission
model, and computes the worst-case voltage and line-flow violations an
attacker can induce by coordinated control of the compromised DERs
(corresponding to Section 5 of the paper).

## Repository layout

```
.
├── Project.toml / Manifest.toml      Julia project environment
│
├── hawaii40.m                        Base 37-bus synthetic Oahu network (MATPOWER format)
├── hawaii40_with_solar.m             Augmented case with compromisable solar inverters added
│
├── add_solar_inverters.jl            Maps exposed inverters to geographically closest 69 kV PV buses
│
├── oahu_solar_inverters.csv          Oahu-localized exposed hosts from the measurement pipeline
├── oahu_solar_inverters.numbers      Spreadsheet copy of the above
├── hosts_lat_lon.csv                 Broader DER host set (ip, lat, lon, label)
│
├── calc_attack_matlab/               Worst-case attack optimization (MATLAB + YALMIP + fmincon)
│   ├── calc_attack.m                     Core optimization (Eq. 3 in the paper)
│   ├── run_calc_attack.m                 Driver: baseline OPF + per-bus / per-line attack loops
│   └── optimal_scurve_k_asym_brute_force.m   AGC sigmoid-fit helper
│
├── calc_attack_julia/                Julia port of the attack optimization (JuMP + Ipopt)
  ├── calc_attack.jl
  ├── run_calc_attack.jl
  └── attack_utils.jl

```

## Pipeline

1. Load the synthetic Oahu 37-bus transmission model (`hawaii40.m`).
2. Augment it with compromisable solar inverters matched from
   `oahu_solar_inverters.csv` to the geographically closest 69 kV PV buses
   (`add_solar_inverters.jl`). The resulting model is `hawaii40_with_solar.m`.
3. Solve a baseline AC OPF to obtain nominal setpoints.
4. For each target bus and line, solve the worst-case attack optimization
   (`calc_attack`), where compromised inverters are free within their
   apparent-power circle and non-compromised generators follow smoothed AGC
   droop and PV/PQ switching.
5. Aggregate the per-target results into the worst-case voltage and
   line-flow distributions reported in the paper.

## Dependencies

**MATLAB path.** MATLAB with the Optimization Toolbox (`fmincon`),
[MATPOWER](https://matpower.org/), and [YALMIP](https://yalmip.github.io/) on
the path.

**Julia path.** Julia 1.9+ with the packages pinned in `Manifest.toml`:
```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running

MATLAB:
```matlab
cd calc_attack_matlab
run_calc_attack
```

Julia:
```
julia --project=. calc_attack_julia/run_calc_attack.jl
```

Both drivers load `hawaii40_with_solar.m`, solve a baseline OPF, and loop over
each bus (voltage targets) and each line (flow targets), and produce the
worst-case histograms.
