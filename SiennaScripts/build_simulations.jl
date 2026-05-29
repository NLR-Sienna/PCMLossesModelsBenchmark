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
    )

    # Create ED model: optimizes dispatch given fixed commitments
    ed_model = make_ptdf_model_without_losses(
        sys_ed;
        device_models = ed_models,
        optimizer = ed_optimizer,
        ptdf = ptdf_ed_used,
        name = "ED",
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
        sys_uc, sys_ed;
        uc_models, ed_models, uc_optimizer, ed_optimizer, ptdf_uc, ptdf_ed
    ) -> Simulation

Build a two-stage UC–ED simulation where:
- UC uses a lossless PTDF model (binary decisions, solved by a MILP optimizer).
- ED uses a PTDF model augmented with quadratic P=I²R losses assuming flat
  voltage (V=1.0 p.u.). Binaries are fixed via SemiContinuousFeedforward from
  UC, so the ED is a pure NLP solvable by Ipopt — no Gurobi required.

Loss = Σₖ Rₖ × (Σⱼ PTDFₖⱼ × Injectionⱼ)²
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
)
    ptdf_uc_used = isnothing(ptdf_uc) ? PTDF(sys_uc) : ptdf_uc
    ptdf_ed_used = isnothing(ptdf_ed) ? PTDF(sys_ed) : ptdf_ed

    sim = build_uc_ed_simulation_with_no_losses(
        sys_uc, sys_ed;
        uc_models, ed_models, uc_optimizer, ed_optimizer,
        ptdf_uc = ptdf_uc_used, ptdf_ed = ptdf_ed_used,
    )

    ed_model = sim.models.decision_models[2]
    update_copperplate_quadratic_loss_approximation_no_voltage!(ed_model, sys_ed, ptdf_ed_used)

    return sim
end