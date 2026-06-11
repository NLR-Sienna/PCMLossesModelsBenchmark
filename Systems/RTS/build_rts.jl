function build_rts_system()
    sys = build_system(PSISystems, "modified_RTS_GMLC_DA_sys"; skip_serialization=true)
    transform_single_time_series!(sys, Hour(1), Hour(1))
    hvdc1 = only(get_components(TwoTerminalGenericHVDCLine, sys))
    set_loss!(hvdc1, LinearCurve(0.0))
    set_reactive_power_limits_from!(hvdc1, (min=-100, max=100))
    set_reactive_power_limits_to!(hvdc1, (min=-100, max=100))
    set_active_power_limits_from!(hvdc1, (min=-100, max=100))
    set_active_power_limits_to!(hvdc1, (min=-100, max=100))
    set_available!(get_component(HydroDispatch, sys, "201_HYDRO_4"), false)
    return sys
end

function set_renewable_costs!(sys, cost_sign::Float64)
    renewables = collect(get_components(RenewableDispatch, sys))
    for (row, gen) in enumerate(renewables)
        new_cost = RenewableGenerationCost(;
            variable=CostCurve(LinearCurve(cost_sign * 0.01 * row)),
        )
        set_operation_cost!(gen, new_cost)
    end
end

function make_uc_template()
    template = ProblemTemplate(NetworkModel(
        PTDFPowerModel;
        use_slacks=true,
        duals=[CopperPlateBalanceConstraint],
        power_flow_evaluation=PowerFlows.ACPowerFlow(; calculate_loss_factors=true, calculate_voltage_stability_factors = true),
    ))
    set_device_model!(template, Line, StaticBranchBounds)
    set_device_model!(template, TapTransformer, StaticBranchBounds)
    set_device_model!(template, ThermalStandard, ThermalBasicUnitCommitment)
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, HydroDispatch, HydroDispatchRunOfRiver)
    set_device_model!(
        template,
        DeviceModel(TwoTerminalGenericHVDCLine, HVDCTwoTerminalDispatch),
    )
    return template
end

"""
    build_rts_model_with_quadratic_losses(sys) -> DecisionModel

Build the UC decision model for the RTS system and augment it with the quadratic
P = I²R loss term, without requiring a prior AC power flow solution.

Steps:
1. Build a PTDF UC model using HiGHS (build! only, no solve)
2. Add PSI-managed loss variable + copper-plate balance update
3. Add quadratic constraint: loss = -Σₖ Rₖ (Σⱼ PTDFₖⱼ Pⱼ)²

Solving the resulting nonconvex MIQP requires Gurobi with NonConvex=2.
"""
function build_rts_model_with_quadratic_losses(sys)
    ptdf = PTDF(sys)
    optimizer = optimizer_with_attributes(
        Gurobi.Optimizer,
        "MIPGap" => 0.001,
        "OutputFlag" => 1,
        "DisplayInterval" => 1,   # print every node
        "LogFile" => "gurobi.log",
        "Presolve" => 2,
        "Heuristics" => 0.3,
        "NonConvex" => 2,
        "MIPFocus" => 2,        # bound improvement
        "Cuts" => 2
    )
    model = PSI.DecisionModel(
        make_uc_template(), sys;
        optimizer=optimizer, name="UC", store_variable_names=true,
    )
    PSI.build!(model, output_dir=mktempdir())
    update_copperplate_quadratic_loss_approximation_no_voltage!(model, sys, ptdf)
    return model
end

"""
    build_rts_uc_models_hv(; voltage_threshold = 100.0) -> Dict

Return UC device models for RTS where `Line` and `TapTransformer` carry a
`filter_function` that excludes branches whose **from-bus** base voltage is at
or below `voltage_threshold` kV from the `NetworkFlowConstraint`.

RTS uses `TapTransformer` (not `Transformer2W` as in CATS).
The full system PTDF must still be passed for correct loss accounting.
"""
function build_rts_uc_models_hv(; voltage_threshold::Float64 = 100.0)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    return Dict(
        Line              => DeviceModel(Line, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        TapTransformer    => DeviceModel(TapTransformer, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicUnitCommitment,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

"""
    build_rts_ed_models_hv(; voltage_threshold = 100.0) -> Dict

Return ED device models for RTS with the same LV branch filter as
`build_rts_uc_models_hv`.
"""
function build_rts_ed_models_hv(; voltage_threshold::Float64 = 100.0)
    filter_fn = x -> PSY.get_base_voltage(PSY.get_from(PSY.get_arc(x))) > voltage_threshold
    return Dict(
        Line              => DeviceModel(Line, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        TapTransformer    => DeviceModel(TapTransformer, StaticBranchBounds;
                                 attributes = Dict("filter_function" => filter_fn)),
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end

"""
    build_rts_ed_models_acopf() -> Dict

Return ED device models for the RTS system using `ACPPowerModel`.
All branches are modeled with `StaticBranchUnbounded` — no `filter_function` —
because `ACPPowerModel` must see every branch for AC feasibility.
RTS uses `TapTransformer` (not `Transformer2W` as in CATS).
"""
function build_rts_ed_models_acopf()
    return Dict(
        Line                       => StaticBranchUnbounded,
        TapTransformer             => StaticBranchUnbounded,
        ThermalStandard            => ThermalBasicDispatch,
        PowerLoad                  => StaticPowerLoad,
        RenewableDispatch          => RenewableFullDispatch,
        HydroDispatch              => HydroDispatchRunOfRiver,
        TwoTerminalGenericHVDCLine => HVDCTwoTerminalLossless,
    )
end