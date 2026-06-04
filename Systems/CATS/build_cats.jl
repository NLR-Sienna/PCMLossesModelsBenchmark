# CATS-specific device models.
# CATS uses Transformer2W (not TapTransformer) and has SynchronousCondenser.
const CATS_UC_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicUnitCommitment,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

const CATS_ED_MODELS = Dict(
    Line => StaticBranchBounds,
    Transformer2W => StaticBranchBounds,
    ThermalStandard => ThermalBasicDispatch,
    PowerLoad => StaticPowerLoad,
    RenewableDispatch => RenewableFullDispatch,
    HydroDispatch => HydroDispatchRunOfRiver,
    SynchronousCondenser => SynchronousCondenserBasicDispatch,
    TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
)

"""
    build_cats_system(cats_json_path) -> System

Load the CATS system from the saved JSON, add two synthetic HVDC links
(Newark–NRS and Metcalf–SanJoseB at ±1000 MW, zero loss), and transform
the time series to 1-hour horizon / 1-hour interval.

The HVDC lines are required to create inter-area shortcuts that the
optimizer can exploit, making circular flows possible.
"""
function build_cats_system(cats_json_path::String)
    sys = System(cats_json_path; runchecks = false)
    transform_single_time_series!(sys, Hour(1), Hour(1))
    add_internal_hvdc!(sys)
    return sys
end

"""
    set_cats_renewable_costs!(sys, cost_sign)

Scale renewable generator costs.
cost_sign > 0 → small positive cost (HVDC loop still profitable → circular flow scenario).
cost_sign < 0 → small negative cost (over-generation incentivized, no HVDC arbitrage → baseline).
"""
function set_cats_renewable_costs!(sys::PSY.System, cost_sign::Float64)
    renewables = collect(PSY.get_components(RenewableDispatch, sys))
    for (row, gen) in enumerate(renewables)
        new_cost = RenewableGenerationCost(;
            variable = CostCurve(LinearCurve(cost_sign * 0.01 * row)),
        )
        PSY.set_operation_cost!(gen, new_cost)
    end
end
