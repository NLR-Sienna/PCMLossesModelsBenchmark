using Statistics
using Printf
using DataFrames

"""
    compare_voltage_stability_factors(sf_a, sf_b; top_n=20) -> DataFrame

Build a per-bus comparison table from two WIDE-format stability factor DataFrames.
`sf_a` = Scenario A (with circular flows), `sf_b` = Scenario B (without).
Both DataFrames must have identical bus columns (same buses in the same column order),
with layout `[DateTime, bus_1, bus_2, ...]`.

Each per-bus value is averaged over all time steps; this handles single-step
simulations trivially and is equally correct for multi-step results.

Returns the top `top_n` rows (default 20) sorted by `|ΔFactor|` descending,
with the following six columns:

- `Bus`                : PSY bus number (integer).
- `Factor_A_CircFlow`  : Mean stability factor for Scenario A.
- `Factor_B_NoCircFlow`: Mean stability factor for Scenario B.
- `Delta_A_minus_B`    : `Factor_A − Factor_B`; positive means the factor is larger
                         in the circular-flow scenario.
- `AbsDelta`           : `|Delta_A_minus_B|`; used for sorting.
- `CircFlows_Smaller`  : `Bool`; `true` when Scenario A's factor is smaller than B's,
                         indicating the bus is closer to a voltage stability limit
                         when circular flows are present.

Note: a *smaller* stability factor means the bus is closer to its voltage stability limit.
"""
function compare_voltage_stability_factors(sf_a::DataFrame, sf_b::DataFrame; top_n::Int = 20)
    bus_cols = names(sf_a)[2:end]
    means_a  = [mean(sf_a[!, c]) for c in bus_cols]
    means_b  = [mean(sf_b[!, c]) for c in bus_cols]

    bus_numbers = parse.(Int, bus_cols)
    delta       = means_a .- means_b
    abs_delta   = abs.(delta)

    tbl = DataFrame(
        Bus                 = bus_numbers,
        Factor_A_CircFlow   = means_a,
        Factor_B_NoCircFlow = means_b,
        Delta_A_minus_B     = delta,
        AbsDelta            = abs_delta,
        CircFlows_Smaller   = means_a .< means_b,
    )
    sort!(tbl, :AbsDelta; rev = true)
    return first(tbl, top_n)
end

"""
    print_stability_comparison(sf_a, sf_b; top_n=20)

Print a formatted per-bus voltage stability factor comparison table between
Scenario A (with circular flows) and Scenario B (without), then summarise
whether circular flows tend to reduce stability factors.

The "✓" flag in the printed output marks buses where Scenario A's stability
factor is strictly smaller than Scenario B's, meaning circular flows are
associated with reduced voltage stability margin at that bus.

The summary line counts how many of the top-`top_n` buses (ranked by `|ΔFactor|`,
default 20) carry the "✓" flag. 
"""
function print_stability_comparison(sf_a::DataFrame, sf_b::DataFrame; top_n::Int = 20)
    comparison = compare_voltage_stability_factors(sf_a, sf_b; top_n = top_n)

    println("\n=== Voltage Stability Factor Comparison (A: circular flows vs B: no circular flows) ===")
    println("  Positive Delta (A−B) → factor LARGER when circular flows present")
    println("  Negative Delta (A−B) → factor SMALLER when circular flows present\n")

    println("  Top $(nrow(comparison)) buses by |ΔFactor| (A−B):")
    println("  $(rpad("Bus", 8)) $(rpad("Factor_A", 12)) $(rpad("Factor_B", 12)) $(rpad("Delta", 12)) CircFlows_Smaller")
    println("  " * "-"^62)
    for row in eachrow(comparison)
        flag = row.CircFlows_Smaller ? "  ✓" : ""
        @printf(
            "  %-8d  %-12.6f  %-12.6f  %+.6f%s\n",
            row.Bus,
            row.Factor_A_CircFlow,
            row.Factor_B_NoCircFlow,
            row.Delta_A_minus_B,
            flag,
        )
    end

    n_smaller = sum(comparison.CircFlows_Smaller)
    n_total   = nrow(comparison)
    total_delta = sum(comparison.Delta_A_minus_B)
    println("\n  Summary: in the top-$n_total buses by |ΔFactor|,")
    println("  $n_smaller / $n_total have a SMALLER stability factor when circular flows are present.")
    println("  Total ΔFactor (A-B) across these buses: $(round(total_delta, sigdigits=4)).")
    if total_delta < 0
        println("  → Circular flows appear to REDUCE voltage stability factors at most affected buses.")
    else
        println("  → Circular flows do NOT consistently reduce stability factors at these buses.")
    end
end
