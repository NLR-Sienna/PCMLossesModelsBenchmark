"""
    build_uc_ed_simulation_with_no_losses(
        sys_uc::PSY.System,
        sys_ed::PSY.System
    ) -> Simulation

Build a two-stage UC-ED simulation without transmission loss modeling.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for Economic Dispatch stage
- `uc_models`: Device models for UC (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC problem (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED problem (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if needed)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if needed)
- `ignore_pf_uc::Bool`: When `true`, the UC stage skips the AC power flow solve
  that normally runs after each optimization step. Recommended for UC because the
  MILP commitment solve is slow and the PF aux variables are rarely needed at the
  UC resolution (default: `false`).
- `ignore_pf_ed::Bool`: When `true`, the ED stage also skips the post-solve AC
  power flow. When `false` (default), an AC power flow is run after each ED solve,
  populating aux variables such as `PowerFlowVoltageStabilityFactors__ACBus`,
  `PowerFlowLossFactors__ACBus`, `PowerFlowVoltageMagnitude__ACBus`, and all branch
  power flows. These aux variables are required for subsequent loss-aware iterations
  (default: `false`).

# Returns
- `Simulation`: Built simulation object ready for execution

# Details
Creates a sequential two-stage simulation structure:

**Stage 1 - Unit Commitment (UC):**
- Determines generator on/off schedules
- Uses PTDF-based network model (DC power flow, lossless)
- Typically runs at hourly resolution
- Solves MILP problem with binary commitment variables

**Stage 2 - Economic Dispatch (ED):**
- Optimizes generation dispatch with fixed commitments from UC
- Also uses PTDF-based network model (lossless)
- Typically runs at sub-hourly resolution (e.g., 5-15 minutes)
- Solves LP or NLP problem (no binary variables)

**Feedforward:**
- SemiContinuousFeedforward links UC to ED
- Passes generator on/off status (OnVariable) from UC to ED
- Ensures ED respects UC commitment decisions
- Affects ActivePowerVariable in ED stage

**Chronology:**
- InterProblemChronology maintains state (e.g., storage levels) between stages
- Ensures temporal consistency across the UC-ED sequence

**Use Case:**
This lossless simulation provides a baseline solution for comparison or
as the initial iteration in iterative loss approximation methods.

# Example
```julia
# Build systems with appropriate time resolutions
sys_uc = build_uc_system("1h")   # Hourly
sys_ed = build_ed_system("15min") # 15-minute

# Create lossless baseline simulation
sim = build_uc_ed_simulation_with_no_losses(sys_uc, sys_ed)
execute!(sim)
results = SimulationResults(sim)
```
"""
function build_uc_ed_simulation_with_no_losses(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
    ignore_pf_uc = false,
    ignore_pf_ed = false,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Create UC model: determines generator commitment schedules
    uc_model = make_ptdf_model_without_losses(
        sys_uc;
        device_models = uc_models,
        optimizer = uc_optimizer,
        ptdf = ptdf_uc_used,
        name = "UC",
        ignore_pf = ignore_pf_uc,
    )

    # Create ED model: optimizes dispatch given fixed commitments
    ed_model = make_ptdf_model_without_losses(
        sys_ed;
        device_models = ed_models,
        optimizer = ed_optimizer,
        ptdf = ptdf_ed_used,
        name = "ED",
        ignore_pf = ignore_pf_ed,
    )

    # Package models into simulation structure
    models = SimulationModels(;
        decision_models = [
            uc_model,
            ed_model,
        ],
    )

    # Define temporal sequence and information flow between stages
    sequence = SimulationSequence(;
        models = models,
        feedforwards = Dict(
            "ED" => [
                # Pass generator on/off status from UC to ED
                # Ensures ED respects UC commitment decisions
                SemiContinuousFeedforward(;
                    component_type = ThermalStandard,
                    source = OnVariable,              # UC commitment decision
                    affected_values = [ActivePowerVariable],  # ED dispatch variable
                ),
            ],
        ),
        # Maintain state variables (e.g., storage) across stages
        ini_cond_chronology = InterProblemChronology(),
    )

    # Create simulation object
    sim = Simulation(;
        name = "Sim",
        steps = 1,                    # Single simulation period
        models = models,
        sequence = sequence,
        simulation_folder = mktempdir(),  # Temporary directory for outputs
    )

    # Build the simulation (construct JuMP models and constraints)
    build!(sim; console_level = Logging.Error)

    return sim
end

"""
    build_uc_ed_simulation_with_uc_linear_ed_quadratic_losses(
        sys_uc::PSY.System,
        sys_ed::PSY.System,
        res_old_uc,
        res_old_ed,
        injection_old_uc,
        injection_old_ed
    ) -> Simulation

Build a two-stage UC-ED simulation with hybrid transmission loss modeling.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for Economic Dispatch stage
- `res_old_uc`: Previous iteration's UC optimization results
- `res_old_ed`: Previous iteration's ED optimization results
- `injection_old_uc`: Previous iteration's UC bus injection values
- `injection_old_ed`: Previous iteration's ED bus injection values
- `uc_models`: Device models for UC (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC problem (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED problem (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if needed)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if needed)

# Returns
- `Simulation`: Built simulation with loss modeling, ready for execution

# Details
Creates a two-stage simulation with differentiated loss modeling strategies:

**UC Stage - Linear Loss Approximation:**
- Loss linearized around previous operating point (injection_old_uc)
- First-order Taylor expansion: Loss ≈ Loss₀ + Σᵢ (∂Loss/∂Pᵢ) × (Pᵢ - Pᵢ,₀)
- Loss factors (∂Loss/∂P) extracted from previous AC power flow
- MILP-compatible formulation (no quadratic terms)
- Allows optimization of binary commitment decisions with loss awareness

**ED Stage - Quadratic Loss Formulation:**
- Accurate P = I²R loss representation
- Loss = Σₖ Rₖ × (Iₖ)² where currents computed via PTDF
- Voltage magnitudes fixed from previous AC power flow
- Requires NLP solver (Ipopt) for quadratic constraints
- More accurate dispatch with fixed commitments from UC

**Rationale for Hybrid Approach:**
1. UC needs MILP solver for commitment decisions → linear losses
2. ED has fixed commitments (no binaries) → can use accurate quadratic losses
3. Balances computational tractability with loss modeling accuracy

**Transmission Constraint Updates:**
Both stages use Fictitious Nodal Demand (FND) method to update branch flow limits:
- Distributes branch losses to endpoint buses as fictitious demands
- Adjusts flow limits via PTDF to reflect loss-induced loading
- Ensures consistency between nodal balance and branch constraints

**Typical Usage Pattern:**
1. First call `build_uc_ed_simulation_with_no_losses` to get initial solution
2. Use those results as input (res_old_*, injection_old_*) to this function
3. Can be called iteratively until convergence

# Example
```julia
# Step 1: Get initial lossless solution
sim0 = build_uc_ed_simulation_with_no_losses(sys_uc, sys_ed)
execute!(sim0)
res0 = SimulationResults(sim0)
res_uc0 = get_decision_problem_results(res0, "UC")
res_ed0 = get_decision_problem_results(res0, "ED")

# Extract injections (requires direct container access)
uc_container = sim0.models.decision_models[1].internal.container
ed_container = sim0.models.decision_models[2].internal.container
optimize!(uc_container.JuMPmodel)
optimize!(ed_container.JuMPmodel)
inj_uc0 = get_injections(uc_container)
inj_ed0 = get_injections(ed_container)

# Step 2: Build simulation with losses
sim = build_uc_ed_simulation_with_uc_linear_ed_quadratic_losses(
    sys_uc, sys_ed, res_uc0, res_ed0, inj_uc0, inj_ed0
)
execute!(sim)
```

# See Also
- `build_uc_ed_simulation_with_no_losses`: Baseline lossless simulation
- `run_iterative_uc_ed_quadratic_loss_simulation`: Iterative wrapper with convergence
"""
function build_uc_ed_simulation_with_uc_linear_ed_quadratic_losses(
    sys_uc::PSY.System,
    sys_ed::PSY.System,
    res_old_uc,
    res_old_ed,
    injection_old_uc,
    injection_old_ed;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Build baseline lossless simulation structure
    # This creates the UC-ED sequence with feedforwards
    sim = build_uc_ed_simulation_with_no_losses(
        sys_uc,
        sys_ed;
        uc_models,
        ed_models,
        uc_optimizer,
        ed_optimizer,
        ptdf_uc = ptdf_uc_used,
        ptdf_ed = ptdf_ed_used,
    )

    # Extract model references for modification
    uc_model = sim.models.decision_models[1]
    ed_model = sim.models.decision_models[2]

    #####################
    ##### UC update #####
    #####################

    # Extract loss parameters from previous UC solution
    loss_factors = get_bus_loss_factors(res_old_uc)      # ∂Loss/∂P at each bus
    total_loss_est = get_total_AC_loss(res_old_uc)       # Total AC losses in MW

    # Add linearized loss approximation to copper plate balance
    # This modifies the copperplate balance: (G - D) + loss_variable = total_loss_est
    # where loss_variable is constrained by:
    #   loss_variable = Σᵢ loss_factors[ᵢ] × (injection_old[ᵢ] - injection[ᵢ])
    # This is a first-order Taylor expansion around the previous operating point
    update_copperplate_loss_approximation!(
        uc_model,
        loss_factors,
        total_loss_est,
        injection_old_uc,
    )

    # Update transmission constraints with fictitious nodal demands
    # This ensures branch flow limits account for the additional loading from losses
    # by distributing losses to buses and propagating via PTDF
    update_transmission_constraints_with_losses!(uc_model, res_old_uc, sys_uc, ptdf_uc_used)

    #####################
    ##### ED update #####
    #####################

    # Add quadratic loss formulation to ED stage
    # Copper plate balance becomes: (G - D) + loss_variable = 0
    # where loss_variable is constrained by:
    #   loss_variable = -Σₖ Rₖ × (Σⱼ PTDFₖⱼ × V_factorsₖⱼ × injectionⱼ)²
    # This represents accurate P = I²R losses with voltage magnitudes from previous flow
    update_copperplate_quadratic_loss_approximation!(
        ed_model,
        sys_ed,
        ptdf_ed_used,
        res_old_ed,
    )

    # Update ED transmission constraints with fictitious nodal demands
    # Same FND method as UC, but based on ED's quadratic loss distribution
    update_transmission_constraints_with_losses!(ed_model, res_old_ed, sys_ed, ptdf_ed_used)

    return sim
end


"""
    build_uc_ed_simulation_with_acopf(
        sys_uc::PSY.System,
        sys_ed::PSY.System
    ) -> Simulation

Build a two-stage UC–ED simulation where the UC stage uses a lossless PTDF
model and the ED stage uses a full AC OPF formulation.

The UC stage must remain PTDF-based (not AC OPF) because binary commitment
variables require a MILP solver, which is incompatible with the nonlinear AC
power flow equations. The ED stage has fixed commitments (no binary variables),
so it can use the full nonlinear AC OPF to capture voltage constraints and
accurate losses.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for the Economic Dispatch stage
- `uc_models`: Device models for UC (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC — must be a MILP solver (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED — must be an NLP solver, e.g. Ipopt
  (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if needed)
- `ptdf_ed`: Pre-computed PTDF matrix for ED; passed to `make_acopf_model` for
  any PTDF-based auxiliary computations (default: nothing, computed if needed)
- `initial_time`: Simulation start timestamp (default: `DateTime("2019-01-01T00:00:00")`)

# Returns
- `Simulation`: Built simulation object ready for execution

# Details
**Stage 1 – Unit Commitment (UC):**
- Lossless PTDF-based network model; binary commitment decisions
- Post-solve AC power flow is skipped (`ignore_pf = true`) — UC is MILP and
  PF aux variables are not needed at the UC resolution

**Stage 2 – Economic Dispatch (ED):**
- Full AC OPF (`make_acopf_model`); voltage magnitudes and angles are
  optimization variables
- Post-solve AC power flow is enabled (`ignore_pf = false`), populating
  aux variables after each ED solve
- Binaries fixed via SemiContinuousFeedforward from UC

# See Also
- `build_uc_ed_simulation_with_acopf_and_uc_linear_losses`: Extends this
  function with a linearized loss approximation in the UC stage
- `build_uc_double_ed_simulation_with_acopf`: Three-stage variant that adds
  an intermediate PTDF-based ED before the AC OPF ED
"""
function build_uc_ed_simulation_with_acopf(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
    initial_time = DateTime("2019-01-01T00:00:00"),
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Create UC model: determines generator commitment schedules
    uc_model = make_ptdf_model_without_losses(
        sys_uc;
        device_models = uc_models,
        optimizer = uc_optimizer,
        ptdf = ptdf_uc_used,
        name = "UC",
        ignore_pf = true,
    )

    # Create ED model: optimizes dispatch given fixed commitments
    ed_model = make_acopf_model(
        sys_ed;
        device_models = ed_models,
        optimizer = ed_optimizer,
        ptdf = ptdf_ed_used,
        name = "ED",
        ignore_pf = false,
    )

    # Package models into simulation structure
    models = SimulationModels(;
        decision_models = [
            uc_model,
            ed_model,
        ],
    )

    # Define temporal sequence and information flow between stages
    sequence = SimulationSequence(;
        models = models,
        feedforwards = Dict(
            "ED" => [
                # Pass generator on/off status from UC to ED
                # Ensures ED respects UC commitment decisions
                SemiContinuousFeedforward(;
                    component_type = ThermalStandard,
                    source = OnVariable,              # UC commitment decision
                    affected_values = [ActivePowerVariable],  # ED dispatch variable
                ),
            ],
        ),
        # Maintain state variables (e.g., storage) across stages
        ini_cond_chronology = InterProblemChronology(),
    )

    # Create simulation object
    sim = Simulation(;
        name = "Sim",
        steps = 1,                    # Single simulation period
        models = models,
        sequence = sequence,
        simulation_folder = mktempdir(),  # Temporary directory for outputs
        initial_time = initial_time,
    )

    # Build the simulation (construct JuMP models and constraints)
    build!(sim; console_level = Logging.Error)

    return sim
end

"""
    build_uc_ed_simulation_with_acopf_and_uc_linear_losses(
        sys_uc::PSY.System,
        sys_ed::PSY.System,
        res_old_uc,
        res_old_ed,
        injection_old_uc
    ) -> Simulation

Build a two-stage UC–ED simulation that combines a linearized loss
approximation in the UC stage with a full AC OPF in the ED stage.

This is the loss-aware variant of `build_uc_ed_simulation_with_acopf`.
The UC stage remains PTDF-based (required for MILP compatibility) but is
augmented with a first-order Taylor expansion of transmission losses around
the previous operating point. The ED stage solves the full AC OPF as before.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for the Economic Dispatch stage
- `res_old_uc`: Previous iteration's UC optimization results (used for loss
  factor extraction — currently unused; loss factors are taken from `res_old_ed`)
- `res_old_ed`: Previous iteration's ED (AC OPF) optimization results; loss
  factors (`∂Loss/∂P`) and total AC losses are extracted from this result
- `injection_old_uc`: Previous iteration's UC bus injection values; used as
  the linearization point for the Taylor-expanded loss expression in UC
- `uc_models`: Device models for UC (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC — must be a MILP solver (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED — must be an NLP solver (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if needed)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if needed)

# Returns
- `Simulation`: Built simulation with UC loss approximation, ready for execution

# Details
Internally calls `build_uc_ed_simulation_with_acopf` to construct the
UC–ED skeleton, then modifies the UC model in-place:

**UC Loss Update:**
- Loss factors and total AC loss are extracted from the previous ED (AC OPF)
  result (`res_old_ed`) via `get_bus_loss_factors` and `get_total_AC_loss_CATS`
- `update_copperplate_loss_approximation!` adds the linearized loss term to
  the UC copper-plate balance constraint

Note: the transmission constraint update
(`update_transmission_constraints_with_losses!`) is currently commented out
for the UC stage; the FND correction is not applied here.

# See Also
- `build_uc_ed_simulation_with_acopf`: Baseline AC OPF simulation (no loss
  approximation in UC)
- `build_uc_ed_simulation_with_uc_linear_ed_quadratic_losses`: Similar hybrid
  approach but with quadratic losses in ED instead of full AC OPF
"""
function build_uc_ed_simulation_with_acopf_and_uc_linear_losses(
    sys_uc::PSY.System,
    sys_ed::PSY.System,
    res_old_uc,
    res_old_ed,
    injection_old_uc;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Build baseline lossless simulation structure
    # This creates the UC-ED sequence with feedforwards
    sim = build_uc_ed_simulation_with_acopf(
        sys_uc,
        sys_ed;
        uc_models,
        ed_models,
        uc_optimizer,
        ed_optimizer,
        ptdf_uc = ptdf_uc_used,
        ptdf_ed = ptdf_ed_used,
    )

    # Extract model references for modification
    uc_model = sim.models.decision_models[1]
    ed_model = sim.models.decision_models[2]

    #####################
    ##### UC update #####
    #####################

    # Extract loss parameters from previous UC solution
    loss_factors = get_bus_loss_factors(res_old_ed; slack_number = "1951")      # ∂Loss/∂P at each bus
    total_loss_est = get_total_AC_loss_CATS(res_old_ed)       # Total AC losses in MW

    # Add linearized loss approximation to copper plate balance
    # This modifies the copperplate balance: (G - D) + loss_variable = total_loss_est
    # where loss_variable is constrained by:
    #   loss_variable = Σᵢ loss_factors[ᵢ] × (injection_old[ᵢ] - injection[ᵢ])
    # This is a first-order Taylor expansion around the previous operating point
    update_copperplate_loss_approximation!(
        uc_model,
        loss_factors,
        total_loss_est,
        injection_old_uc,
    )

    # Update transmission constraints with fictitious nodal demands
    # This ensures branch flow limits account for the additional loading from losses
    # by distributing losses to buses and propagating via PTDF
    # update_transmission_constraints_with_losses!(uc_model, res_old_uc, sys_uc, ptdf_uc_used)

    return sim
end


"""
    build_uc_double_ed_simulation_with_acopf(
        sys_uc::PSY.System,
        sys_ed::PSY.System
    ) -> Simulation

Build a three-stage simulation: UC (MILP, lossless PTDF) → EDPTDF (LP,
lossless PTDF) → ED (NLP, full AC OPF).

The intermediate PTDF-based ED stage ("EDPTDF") sits between the binary
commitment UC and the nonlinear AC OPF ED. It serves as a warm-start or
feasibility-screening step before the full AC solve: commitment decisions
propagate from UC to both EDPTDF and ED via separate SemiContinuousFeedforward
links.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System shared by EDPTDF and ED stages
- `uc_models`: Device models for UC (default: DEFAULT_UC_MODELS)
- `ed_ptdf_models`: Device models for the intermediate PTDF ED (default: DEFAULT_ED_MODELS)
- `ed_models`: Device models for the AC OPF ED (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC — must be a MILP solver (default: DEFAULT_MILP_OPTIMIZER)
- `ed_ptdf_optimizer`: Optimizer for the PTDF ED stage (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for the AC OPF ED — must be an NLP solver
  (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if needed)
- `ptdf_ed`: Pre-computed PTDF matrix for EDPTDF and ED (default: nothing, computed if needed)

# Returns
- `Simulation`: Built three-stage simulation object ready for execution

# Details
**Stage 1 – Unit Commitment (UC):**
- Lossless PTDF model; binary commitment decisions
- Post-solve AC power flow skipped (`ignore_pf = true`)

**Stage 2 – EDPTDF (intermediate LP Economic Dispatch):**
- Lossless PTDF model; commitments fixed by feedforward from UC
- Post-solve AC power flow skipped (`ignore_pf = true`)

**Stage 3 – ED (AC OPF Economic Dispatch):**
- Full AC OPF; voltage magnitudes and angles are optimization variables
- Post-solve AC power flow enabled (`ignore_pf = false`)
- Commitments fixed by a separate feedforward directly from UC

Both EDPTDF and ED receive the `SemiContinuousFeedforward` (OnVariable →
ActivePowerVariable) independently from UC.

# See Also
- `build_uc_ed_simulation_with_acopf`: Two-stage variant (no intermediate PTDF ED)
- `build_uc_double_ed_simulation_with_acopf_and_uc_linear_ed_quadratic_losses`:
  Extends this structure with loss modeling in UC and ED stages
"""
function build_uc_double_ed_simulation_with_acopf(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_ptdf_models = DEFAULT_ED_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_ptdf_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Create UC model: determines generator commitment schedules
    uc_model = make_ptdf_model_without_losses(
        sys_uc;
        device_models = uc_models,
        optimizer = uc_optimizer,
        ptdf = ptdf_uc_used,
        name = "UC",
        ignore_pf = true,
    )

    ed_ptdf_model = make_ptdf_model_without_losses(
        sys_ed;
        device_models = ed_ptdf_models,
        optimizer = ed_ptdf_optimizer,
        ptdf = ptdf_ed_used,
        name = "EDPTDF",
        ignore_pf = true,
    )

    # Create ED model: optimizes dispatch given fixed commitments
    ed_model = make_acopf_model(
        sys_ed;
        device_models = ed_models,
        optimizer = ed_optimizer,
        ptdf = ptdf_ed_used,
        name = "ED",
        ignore_pf = false,
    )

    # Package models into simulation structure
    models = SimulationModels(;
        decision_models = [
            uc_model,
            ed_ptdf_model,
            ed_model,
        ],
    )

    # Define temporal sequence and information flow between stages
    sequence = SimulationSequence(;
        models = models,
        feedforwards = Dict(
            "EDPTDF" => [
                # Pass generator on/off status from UC to EDPTDF
                # Ensures EDPTDF respects UC commitment decisions
                SemiContinuousFeedforward(;
                    component_type = ThermalStandard,
                    source = OnVariable,              # UC commitment decision
                    affected_values = [ActivePowerVariable],  # ED dispatch variable
                ),
            ],
            "ED" => [
                # Pass generator on/off status from UC to ED
                # Ensures ED respects UC commitment decisions
                SemiContinuousFeedforward(;
                    component_type = ThermalStandard,
                    source = OnVariable,              # UC commitment decision
                    affected_values = [ActivePowerVariable],  # ED dispatch variable
                ),
            ],
        ),
        # Maintain state variables (e.g., storage) across stages
        ini_cond_chronology = InterProblemChronology(),
    )

    # Create simulation object
    sim = Simulation(;
        name = "Sim",
        steps = 1,                    # Single simulation period
        models = models,
        sequence = sequence,
        simulation_folder = mktempdir(),  # Temporary directory for outputs
    )

    # Build the simulation (construct JuMP models and constraints)
    build!(sim; console_level = Logging.Error)

    return sim
end

"""
    build_uc_double_ed_simulation_with_acopf_and_uc_linear_ed_quadratic_losses(
        sys_uc::PSY.System,
        sys_ed::PSY.System,
        res_old_uc,
        res_old_ed,
        injection_old_uc
    ) -> Simulation

Build a UC–AC OPF ED simulation augmented with a linearized loss
approximation in the UC stage.

**Note:** Despite the "double_ed" in the name, this function currently
delegates to `build_uc_ed_simulation_with_acopf` (two-stage skeleton),
not the three-stage `build_uc_double_ed_simulation_with_acopf`. The ED
quadratic loss update section is present as a placeholder but not yet
populated. The function is effectively equivalent to
`build_uc_ed_simulation_with_acopf_and_uc_linear_losses` for the UC update
path.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for the Economic Dispatch stage
- `res_old_uc`: Previous iteration's UC results (signature compatibility;
  loss factors are currently sourced from `res_old_ed`)
- `res_old_ed`: Previous iteration's ED (AC OPF) results; provides loss
  factors and total AC loss estimate for the UC linearization
- `injection_old_uc`: Previous iteration's UC bus injections; used as
  the linearization point for the loss Taylor expansion
- `uc_models`: Device models for UC (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for UC — must be a MILP solver (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for ED — must be an NLP solver (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if needed)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if needed)

# Returns
- `Simulation`: Built simulation with UC loss approximation, ready for execution

# See Also
- `build_uc_ed_simulation_with_acopf_and_uc_linear_losses`: Functionally
  equivalent two-stage version
- `build_uc_double_ed_simulation_with_acopf`: Three-stage baseline without
  loss modeling
"""
function build_uc_double_ed_simulation_with_acopf_and_uc_linear_ed_quadratic_losses(
    sys_uc::PSY.System,
    sys_ed::PSY.System,
    res_old_uc,
    res_old_ed,
    injection_old_uc;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
)
    # Use provided PTDF matrices or compute them if not provided
    if isnothing(ptdf_uc)
        ptdf_uc_used = PTDF(sys_uc)
    else
        ptdf_uc_used = ptdf_uc
    end
    if isnothing(ptdf_ed)
        ptdf_ed_used = PTDF(sys_ed)
    else
        ptdf_ed_used = ptdf_ed
    end

    # Build baseline lossless simulation structure
    # This creates the UC-ED sequence with feedforwards
    sim = build_uc_ed_simulation_with_acopf(
        sys_uc,
        sys_ed;
        uc_models,
        ed_models,
        uc_optimizer,
        ed_optimizer,
        ptdf_uc = ptdf_uc_used,
        ptdf_ed = ptdf_ed_used,
    )

    # Extract model references for modification
    uc_model = sim.models.decision_models[1]
    ed_model = sim.models.decision_models[2]

    #####################
    ##### UC update #####
    #####################

    # Extract loss parameters from previous UC solution
    loss_factors = get_bus_loss_factors(res_old_ed; slack_number = "1951")      # ∂Loss/∂P at each bus
    total_loss_est = get_total_AC_loss_CATS(res_old_ed)       # Total AC losses in MW

    # Add linearized loss approximation to copper plate balance
    # This modifies the copperplate balance: (G - D) + loss_variable = total_loss_est
    # where loss_variable is constrained by:
    #   loss_variable = Σᵢ loss_factors[ᵢ] × (injection_old[ᵢ] - injection[ᵢ])
    # This is a first-order Taylor expansion around the previous operating point
    update_copperplate_loss_approximation!(
        uc_model,
        loss_factors,
        total_loss_est,
        injection_old_uc,
    )

    # Update transmission constraints with fictitious nodal demands
    # This ensures branch flow limits account for the additional loading from losses
    # by distributing losses to buses and propagating via PTDF
    # update_transmission_constraints_with_losses!(uc_model, res_old_uc, sys_uc, ptdf_uc_used)

    #####################
    ##### ED update #####
    #####################

    return sim
end

"""
    build_uc_ed_simulation_with_ed_quadratic_losses_no_voltage(
        sys_uc::PSY.System,
        sys_ed::PSY.System
    ) -> Simulation

Build a two-stage UC–ED simulation where the ED stage uses quadratic P=I²R
losses computed under a flat-voltage (V=1.0 p.u.) assumption.

The "no_voltage" in the name refers specifically to the **loss formulation**:
voltage magnitudes are not corrected during the quadratic loss computation,
so the loss expression simplifies to

    Loss = Σₖ Rₖ × (Σⱼ PTDFₖⱼ × Injectionⱼ)²

This is distinct from the **post-solve AC power flow** controlled by
`ignore_pf_ed`. Even with "no voltage in losses", the ED stage can still
run a full AC power flow after each solve (when `ignore_pf_ed = false`) to
populate aux variables for monitoring and downstream iterations. The flat-
voltage assumption only affects how losses enter the copper-plate balance
during optimization, not the post-solve diagnostics.

# Arguments
- `sys_uc::PSY.System`: PowerSystems.jl System for the Unit Commitment stage
- `sys_ed::PSY.System`: PowerSystems.jl System for the Economic Dispatch stage
- `uc_models`: Device models for UC (default: DEFAULT_UC_MODELS)
- `ed_models`: Device models for ED (default: DEFAULT_ED_MODELS)
- `uc_optimizer`: Optimizer for the UC problem (default: DEFAULT_MILP_OPTIMIZER)
- `ed_optimizer`: Optimizer for the ED problem (default: DEFAULT_NLP_OPTIMIZER)
- `ptdf_uc`: Pre-computed PTDF matrix for UC (default: nothing, computed if needed)
- `ptdf_ed`: Pre-computed PTDF matrix for ED (default: nothing, computed if needed)
- `ignore_pf_uc::Bool`: When `true`, the UC stage skips the post-solve AC power
  flow. Recommended for UC because the MILP solve is slow and PF aux variables
  are rarely needed at UC resolution (default: `false`).
- `ignore_pf_ed::Bool`: When `true`, the ED stage skips the post-solve AC power
  flow. When `false` (default), an AC power flow runs after each ED solve,
  populating aux variables including `PowerFlowVoltageStabilityFactors__ACBus`,
  `PowerFlowLossFactors__ACBus`, `PowerFlowVoltageMagnitude__ACBus`, and branch
  power flows. These aux variables are needed for downstream loss-aware iterations
  even though they do not feed back into the flat-voltage loss expression used
  during the ED solve itself.

# Returns
- `Simulation`: Built simulation object ready for execution

# Details
**Stage 1 – Unit Commitment (UC):**
- Lossless PTDF-based network model; binary commitment decisions
- Solved by a MILP optimizer (e.g., Gurobi, HiGHS)

**Stage 2 – Economic Dispatch (ED):**
- PTDF-based network model augmented with quadratic losses (flat V=1 p.u.)
- Binaries fixed via SemiContinuousFeedforward from UC → pure NLP, no MILP solver required
- Uses `update_copperplate_quadratic_loss_approximation_no_voltage_untracked!`
  (untracked variant because HDF5 result store has no pre-allocated slot for
  variables added after build time)

# See Also
- `build_uc_ed_simulation_with_no_losses`: Baseline lossless simulation (also used
  internally to construct the UC–ED skeleton)
- `build_uc_ed_simulation_with_uc_linear_ed_quadratic_losses`: Hybrid simulation
  where UC uses a linearized loss approximation and ED uses voltage-corrected
  quadratic losses
"""
function build_uc_ed_simulation_with_ed_quadratic_losses_no_voltage(
    sys_uc::PSY.System,
    sys_ed::PSY.System;
    uc_models = DEFAULT_UC_MODELS,
    ed_models = DEFAULT_ED_MODELS,
    uc_optimizer = DEFAULT_MILP_OPTIMIZER,
    ed_optimizer = DEFAULT_NLP_OPTIMIZER,
    ptdf_uc = nothing,
    ptdf_ed = nothing,
    ignore_pf_uc = false,
    ignore_pf_ed = false,
)
    ptdf_uc_used = isnothing(ptdf_uc) ? PTDF(sys_uc) : ptdf_uc
    ptdf_ed_used = isnothing(ptdf_ed) ? PTDF(sys_ed) : ptdf_ed

    sim = build_uc_ed_simulation_with_no_losses(
        sys_uc, sys_ed;
        uc_models, ed_models, uc_optimizer, ed_optimizer,
        ptdf_uc = ptdf_uc_used, ptdf_ed = ptdf_ed_used,
        ignore_pf_uc = ignore_pf_uc,
        ignore_pf_ed  = ignore_pf_ed,
    )

    ed_model = sim.models.decision_models[2]
    # Use the untracked variant: HDF5 store has no slot for post-build variables.
    update_copperplate_quadratic_loss_approximation_no_voltage_untracked!(ed_model, sys_ed, ptdf_ed_used)

    return sim
end