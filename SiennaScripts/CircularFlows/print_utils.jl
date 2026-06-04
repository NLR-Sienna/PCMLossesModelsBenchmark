using Statistics
using Printf
using DataFrames

"""
    compare_voltage_stability_factors(sf_a, sf_b; top_n=20) -> DataFrame

Build a per-bus comparison table from two WIDE-format stability factor DataFrames.
sf_a = Scenario A (with circular flows), sf_b = Scenario B (without).
Columns: [DateTime, bus_1, bus_2, ...]. Returns a DataFrame sorted by |Δ| descending.
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
    println("\n  Summary: in the top-$n_total buses by |ΔFactor|,")
    println("  $n_smaller / $n_total have a SMALLER stability factor when circular flows are present.")
    if n_smaller > n_total ÷ 2
        println("  → Circular flows appear to REDUCE voltage stability factors at most affected buses.")
    else
        println("  → Circular flows do NOT consistently reduce stability factors at these buses.")
    end
end
