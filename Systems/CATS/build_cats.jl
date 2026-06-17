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
    build_cats_uc_models_hv(; voltage_threshold = 100.0) -> Dict

Return UC device models for CATS where `Line` and `Transformer2W` carry a
`filter_function` that excludes branches whose **from-bus** base voltage is at
or below `voltage_threshold` kV from the `NetworkFlowConstraint`.

The full system PTDF must still be passed to the PSI `NetworkModel`; this filter
only removes the per-branch flow-limit constraints for LV branches — it does NOT
remove their contribution to the quadratic loss term.
"""
function build_cats_uc_models_hv(; voltage_threshold::Float64 = 100.0, bounded = false)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    if bounded
        branch_model = StaticBranchBounds
    else
        branch_model = StaticBranchUnbounded
    end
    return Dict(
        Line          => DeviceModel(Line, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        Transformer2W => DeviceModel(Transformer2W, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicUnitCommitment,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

"""
    build_cats_ed_models_hv(; voltage_threshold = 100.0) -> Dict

Return ED device models for CATS with the same LV branch filter as
`build_cats_uc_models_hv`. Includes `SynchronousCondenser` which is
present only in the CATS ED model.
"""
function build_cats_ed_models_hv(; voltage_threshold::Float64 = 100.0, bounded = false)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    if bounded
        branch_model = StaticBranchBounds
    else
        branch_model = StaticBranchUnbounded
    end
    return Dict(
        Line          => DeviceModel(Line, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        Transformer2W => DeviceModel(Transformer2W, branch_model;
                             attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        SynchronousCondenser       => SynchronousCondenserBasicDispatch,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

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

function scale_cats_loads!(sys::PSY.System, scale_factor::Float64)
    loads = collect(PSY.get_components(PowerLoad, sys))
    for load in loads
        max_p = PSY.get_max_active_power(load)
        p = PSY.get_active_power(load)
        PSY.set_active_power!(load, p * scale_factor)
        PSY.set_max_active_power!(load, max_p * scale_factor)
    end
end

"""
    build_cats_ed_models_acopf() -> Dict

Return ED device models for the CATS system using `ACPPowerModel`.
All branches are modeled with `StaticBranchUnbounded` — no `filter_function` —
because `ACPPowerModel` must see every branch for AC feasibility.
Losses are implicit in the nonlinear AC formulation; no separate loss term is needed.
"""
function build_cats_ed_models_acopf(;bounded = false)
    if bounded
        branch_model = StaticBranchBounds
    else
        branch_model = StaticBranchUnbounded
    end
    return Dict(
        Line                       => branch_model,
        Transformer2W              => branch_model,
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        SynchronousCondenser       => SynchronousCondenserBasicDispatch,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end
