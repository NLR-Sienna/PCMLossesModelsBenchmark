# needed to add this due to quirks of my local Xpress installation.
ENV["XPRESSDIR"] = "/Users/lkiernan/Documents/Xpress"

# right now SynchronousCondenserBasicDispatch isn't implement on main. Remove Pkg.add line
# once https://github.com/NREL-Sienna/PowerSimulations.jl/pull/1507 is merged.
using Pkg
Pkg.add(url="https://github.com/NREL-Sienna/PowerSimulations.jl.git", rev="rh/add_syncon_model")
using PowerSimulations

using HydroPowerSimulations
using Ipopt
using Xpress
using Dates
using JuMP
using PowerFlows
using PowerNetworkMatrices

const PSI = PowerSimulations


# git clone https://github.com/NREL-Sienna/CATS-CaliforniaTestSystem/tree/lk/redo-Sienna-scripts-rebase
# make sure you're on the branch lk/redo-Sienna-scripts-rebase
CATS_DIR = "/Users/lkiernan/Documents/julia/CATS-project-2/CATS-CaliforniaTestSystem/Sienna/"
include(joinpath(CATS_DIR, "build_CATS.jl"))
# with quadratic cost functions, first time step's optimization problem takes 30+ minutes.
# worth looking into replacing quadratic costs with piecewise linear approximations.
system = build_CATS_system(first_order = true)


transform_single_time_series!(
    system,
    Hour(1),  # horizon
    Hour(1),  # interval
);
const PNM = PowerNetworkMatrices

# PSI has a hard-coded assumption that parallel lines have identical impedances.
# will open PR to remove that assumption.
arc_to_line_impedance = Dict{Tuple{Int, Int}, ComplexF64}()
for line in get_components(Line, system)
    arc = get_arc(line)
    arc_tuple = (arc.from.number, arc.to.number)
    if haskey(arc_to_line_impedance, arc_tuple)
        set_r!(line, real(arc_to_line_impedance[arc_tuple]))
        set_x!(line, imag(arc_to_line_impedance[arc_tuple]))
    else
        @assert !haskey(arc_to_line_impedance, reverse(arc_tuple))
    end
    arc_to_line_impedance[arc_tuple] = get_r(line) + im * get_x(line)
end

# UC PTDF
template_uc = ProblemTemplate(
    NetworkModel(PTDFPowerModel;
        use_slacks=true,
        duals=[CopperPlateBalanceConstraint]
    )
)

set_device_model!(template_uc, ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
set_device_model!(template_uc, Line, StaticBranch)
set_device_model!(template_uc, Transformer2W, StaticBranch)
solver_xpress = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.01, "MAXTIME" => 60*60)

problem_uc = DecisionModel(
    template_uc,
    system;
    optimizer=solver_xpress,
    optimizer_solve_log_print=true,
    calculate_conflict=true,
    name="UC"
)

# ED ACOPF
network_model_ed = NetworkModel(ACPPowerModel; use_slacks=true, power_flow_evaluation = PowerFlows.ACPowerFlow(; calculate_loss_factors=true, calculate_voltage_stability_factors = true))
template_ed = ProblemTemplate(network_model_ed)
set_device_model!(template_ed, ThermalStandard, ThermalBasicDispatch)
set_device_model!(template_ed, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_ed, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_ed, PowerLoad, StaticPowerLoad)
set_device_model!(template_ed, Line, StaticBranchUnbounded)
set_device_model!(template_ed, Transformer2W, StaticBranchUnbounded)
set_device_model!(template_ed, SynchronousCondenser, SynchronousCondenserBasicDispatch)

solver_ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)
problem_ed = DecisionModel(
    template_ed, 
    system; 
    optimizer = solver_ipopt,
    optimizer_solve_log_print = true, 
    name = "ED"
)

models = SimulationModels(;
    decision_models = [
        problem_uc,
        problem_ed,
    ],
)

sequence = SimulationSequence(;
    models = models,
    feedforwards = Dict(
        "ED" => [
            SemiContinuousFeedforward(;
                component_type = ThermalStandard,
                source = OnVariable,
                affected_values = [ActivePowerVariable],
            ),
        ],
    ),
    ini_cond_chronology = InterProblemChronology(),
)

sim = Simulation(;
    name = "no_cache",
    steps = 5, # adjust appropriately. Runtime is 1-2 minutes per time step.
    models = models,
    sequence = sequence,
    simulation_folder = mktempdir(),
)


build!(sim)
execute!(sim)

results_sim = SimulationResults(sim)
results_ed = get_decision_problem_results(results_sim, "ED")

# first time step fails to converge, so showing results from second time step
read_aux_variable(results_ed, "PowerFlowLossFactors__ACBus")[DateTime("2019-01-01T01:00:00")]
read_aux_variable(results_ed,  "PowerFlowVoltageStabilityFactors__ACBus")[DateTime("2019-01-01T01:00:00")]

