###########################################################################################
# Auxiliary methods for processing loss factors terms #####################################
# Methods attempt to return always vectors or matrices rather than DataFrames #############
# The axes are returned in a similar fashion on how DenseAxisArrays are used in container #
# That is the first axes is the name/number/component and second axis is the timestep #####
###########################################################################################

"""
    get_bus_ax(res) -> Vector{Int}

Extract bus numbers from optimization results as an integer vector.

# Arguments
- `res`: SimulationResults or DecisionModelResults containing power balance expressions

# Returns
- Vector of bus numbers (as integers) corresponding to the system buses

# Details
Reads the ActivePowerBalance expression for AC buses and extracts bus identifiers
from the column names, converting them from strings to integers.
"""
function get_bus_ax(res)
    injection = read_expression(res, "ActivePowerBalance__ACBus"; table_format = TableFormat.WIDE)
    bus_numbers_string = names(injection)[2:end]  # Skip first column (timestamp)
    return parse.(Int, bus_numbers_string)
end

"""
    get_bus_injection(res) -> Matrix{Float64}

Retrieve net active power injection at each bus across all time periods.

# Arguments
- `res`: SimulationResults or DecisionModelResults

# Returns
- Matrix of size (num_buses × num_timesteps) containing net injections in MW
  - Row i: injections for bus i
  - Column t: injections at time t

# Details
Net injection = Generation - Demand at each bus.
Positive values indicate net generation, negative values indicate net load.
"""
function get_bus_injection(res)
    injection = read_expression(res, "ActivePowerBalance__ACBus"; table_format = TableFormat.WIDE)
    # Transpose to get (buses × time) format
    return copy(transpose(Matrix{Float64}(injection[!, 2:end])))
end

"""
    get_bus_loss_factors(res; slack_number = "113") -> Matrix{Float64}

Compute bus-specific loss sensitivity factors (∂Loss/∂Injection) from power flow results.

# Arguments
- `res`: SimulationResults or DecisionModelResults containing loss factors
- `slack_number`: Bus number of the slack/reference bus (default: "113")

# Returns
- Matrix of size (num_buses × num_timesteps) containing loss factors
  - loss_factors[i,t] = sensitivity of total system losses to injection at bus i, time t
  - Slack bus factors are set to 0.0 (reference)

# Details
Loss factors represent how much total system losses change per unit change in injection.
They are computed as 1.0 + delivery_factors, where delivery factors account for
the marginal transmission losses. The slack bus is excluded from loss calculations.
"""
function get_bus_loss_factors(res; slack_number = "113")
    delivery_factors = read_aux_variable(res, "PowerFlowLossFactors__ACBus"; table_format = TableFormat.WIDE)
    bus_numbers_string = names(delivery_factors)[2:end]
    ix_slack = findfirst(bus_numbers_string.== slack_number)
    
    # Convert delivery factors to loss factors: LF = 1 + DF
    loss_factors =  1.0 .+ Matrix{Float64}(delivery_factors[!, 2:end])
    
    # Set slack bus loss factor to zero (reference)
    loss_factors[:, ix_slack] .= 0.0
    
    return copy(transpose(loss_factors))
end

function get_bus_loss_factors(res::PSI.SimulationProblemResults{PowerSimulations.DecisionModelSimulationResults}; slack_number = "113")
    delivery_factors = read_realized_aux_variable(res, "PowerFlowLossFactors__ACBus"; table_format = TableFormat.WIDE)
    bus_numbers_string = names(delivery_factors)[2:end]
    ix_slack = findfirst(bus_numbers_string.== slack_number)
    
    # Convert delivery factors to loss factors: LF = 1 + DF
    loss_factors =  1.0 .+ Matrix{Float64}(delivery_factors[!, 2:end])
    
    # Set slack bus loss factor to zero (reference)
    loss_factors[:, ix_slack] .= 0.0
    
    return copy(transpose(loss_factors))
end

"""
    get_total_R_per_arc(sys, bus_no_from, bus_no_to) -> Float64

Compute total series resistance for all parallel branches connecting two buses.

# Arguments
- `sys`: PowerSystems.jl System object
- `bus_no_from`: Source bus number
- `bus_no_to`: Destination bus number

# Returns
- Total resistance in per-unit (sum of all parallel branches)
- Returns 0.0 if no branches found

# Details
For parallel branches between the same buses, this sums the individual resistances.
Note: This is the total resistance, not the equivalent parallel resistance.
Only includes available (in-service) branches.
"""
function get_total_R_per_arc(sys, bus_no_from, bus_no_to)
    # Find all available branches connecting the two buses
    all_branches = get_components(x -> get_available(x) && (x.arc.from.number == bus_no_from) && (x.arc.to.number == bus_no_to), ACBranch, sys)
    if isempty(all_branches)
        println("No branches found from bus $bus_no_from to bus $bus_no_to")
        return 0.0
    end
    # Sum resistances of all parallel branches
    return sum(get_r(branch) for branch in all_branches)
end

"""
    get_total_X_per_arc(sys, bus_no_from, bus_no_to) -> Float64

Compute total series reactance for all parallel branches connecting two buses.

# Arguments
- `sys`: PowerSystems.jl System object
- `bus_no_from`: Source bus number
- `bus_no_to`: Destination bus number

# Returns
- Total reactance in per-unit (sum of all parallel branches)
- Returns 0.0 if no branches found

# Details
For parallel branches between the same buses, this sums the individual reactances.
Note: This is the total reactance, not the equivalent parallel reactance.
Only includes available (in-service) branches.
"""
function get_total_X_per_arc(sys, bus_no_from, bus_no_to)
    # Find all available branches connecting the two buses
    all_branches = get_components(x -> get_available(x) && (x.arc.from.number == bus_no_from) && (x.arc.to.number == bus_no_to), ACBranch, sys)
    if isempty(all_branches)
        println("No branches found from bus $bus_no_from to bus $bus_no_to")
        return 0.0
    end
    # Sum reactances of all parallel branches
    return sum(get_x(branch) for branch in all_branches)
end

"""
    get_power_flow_voltage_mag(res) -> Matrix{Float64}

Retrieve voltage magnitudes at all buses from power flow results.

# Arguments
- `res`: SimulationResults or DecisionModelResults with power flow variables

# Returns
- Matrix of size (num_buses × num_timesteps) containing voltage magnitudes in per-unit
  - Row i: voltage magnitude for bus i
  - Column t: voltages at time t

# Details
Voltage magnitudes are computed from AC power flow analysis and are typically
close to 1.0 p.u. for well-operated systems.
"""
function get_power_flow_voltage_mag(res)
    pf_bus = read_aux_variable(res, "PowerFlowVoltageMagnitude__ACBus"; table_format = TableFormat.WIDE)
    # Transpose to get (buses × time) format
    return copy(transpose(Matrix{Float64}(pf_bus[!, 2:end])))
end

function get_power_flow_voltage_mag(res::PSI.SimulationProblemResults{PowerSimulations.DecisionModelSimulationResults})
    pf_bus = read_realized_aux_variable(res, "PowerFlowVoltageMagnitude__ACBus"; table_format = TableFormat.WIDE)
    # Transpose to get (buses × time) format
    return copy(transpose(Matrix{Float64}(pf_bus[!, 2:end])))
end

"""
    get_power_flow_arc_voltage_mag(res, ptdf) -> Matrix{Float64}

Compute representative voltage magnitude for each transmission arc.

# Arguments
- `res`: SimulationResults or DecisionModelResults
- `ptdf`: PTDF matrix object containing arc definitions

# Returns
- Matrix of size (num_arcs × num_timesteps) containing arc voltages in per-unit
  - Row k: representative voltage for arc k
  - Column t: voltages at time t

# Details
For each arc (branch), uses the maximum voltage magnitude between the two
connected buses. This conservative approach ensures loss calculations
account for the higher voltage conditions.
"""
function get_power_flow_arc_voltage_mag(res, ptdf)
    bus_lookup = ptdf.lookup[1]
    arc_axes = axes(ptdf, 2)  # Get (from_bus, to_bus) pairs
    V_bus = get_power_flow_voltage_mag(res)
    V_line = zeros(length(arc_axes), size(V_bus, 2))
    
    # Use maximum voltage of the two buses for each arc
    for (ix, (bus_no_from, bus_no_to)) in enumerate(arc_axes)
        V_line[ix, :] = max(V_bus[bus_lookup[bus_no_from], :], V_bus[bus_lookup[bus_no_to], :])
    end
    return V_line
end

"""
    get_RX_vector(sys, ptdf) -> Tuple{Vector{Float64}, Vector{Float64}}

Extract resistance and reactance vectors for all transmission arcs.

# Arguments
- `sys`: PowerSystems.jl System object
- `ptdf`: PTDF matrix object defining arc topology

# Returns
- Tuple (R, X) where:
  - R: Vector of resistances (in per-unit) for each arc
  - X: Vector of reactances (in per-unit) for each arc

# Details
Arcs represent unique (from_bus, to_bus) pairs from the PTDF matrix.
For each arc, computes the total R and X including all parallel branches.
"""
function get_RX_vector(sys, ptdf)
    AA = axes(ptdf, 2)  # Arc Axes: Line + Transformer Arcs (unique pairs)
    total_br = length(AA)
    R = zeros(total_br)
    X = zeros(total_br)
    
    # Compute R and X for each arc
    for (ix, (bus_no_from, bus_no_to)) in enumerate(AA)
        R[ix] = get_total_R_per_arc(sys, bus_no_from, bus_no_to)
        X[ix] = get_total_X_per_arc(sys, bus_no_from, bus_no_to)
    end
    return R, X
end

"""
    get_AC_dLoss_dP(res, sys) -> Matrix{Float64}

Compute AC loss sensitivity factors (∂Loss/∂P) accounting for voltage variations.

# Arguments
- `res`: SimulationResults or DecisionModelResults with power flow data
- `sys`: PowerSystems.jl System object

# Returns
- Matrix of size (num_buses × num_timesteps) containing AC loss factors
  - loss_p[j,t] = ∂(Total Loss)/∂(Injection at bus j) at time t

# Details
Computes the gradient of total system active power losses with respect to
bus injections using AC power flow formulation:

    ∂Loss/∂Pⱼ = Σₖ 2·Rₖ·(Vₖ/Vⱼ)·PTDFₖⱼ·Σᵢ (Vₖ/Vᵢ)·PTDFₖᵢ·Pᵢ

where:
- Rₖ: resistance of arc k
- Vⱼ, Vᵢ: voltage magnitudes at buses
- Vₖ: representative voltage for arc k
- PTDFₖⱼ: power transfer distribution factor
- Pᵢ: injection at bus i

This formulation accounts for voltage magnitude effects on losses,
making it more accurate than DC approximation for systems with significant
voltage variations.
"""
function get_AC_dLoss_dP(res, sys)
    ptdf = PTDF(sys)
    V_bus = get_power_flow_voltage_mag(res)
    V_line = get_power_flow_arc_voltage_mag(res, ptdf)
    T_length = size(V_bus, 2)
    bus_length = length(axes(ptdf, 1))
    arcs_length = length(axes(ptdf, 2))
    loss_p = zeros(bus_length, T_length)
    injection = get_bus_injection(res)
    R, _ = get_RX_vector(sys, ptdf)
    
    # Compute AC loss sensitivity for each bus and time period
    # This accounts for voltage-dependent losses
    for j in 1:bus_length
        for t in 1:T_length
            loss_p[j, t] = sum(2 * R[k] * V_line[k] / V_bus[j, t] * ptdf[k, j] * (sum(V_line[k, t] / V_bus[i, t] * ptdf[k, i] * injection[i, t] for i in 1:bus_length)) for k in 1:arcs_length)
        end
    end
    return loss_p
end

"""
    get_DC_dLoss_dP(res, sys) -> Matrix{Float64}

Compute DC loss sensitivity factors (∂Loss/∂P) using constant voltage assumption.

# Arguments
- `res`: SimulationResults or DecisionModelResults
- `sys`: PowerSystems.jl System object

# Returns
- Matrix of size (num_buses × num_timesteps) containing DC loss factors
  - loss_p[i,t] = ∂(Total Loss)/∂(Injection at bus i) at time t

# Details
Computes the gradient of total system active power losses with respect to
bus injections using DC power flow approximation:

    ∂Loss/∂Pᵢ = Σₖ 2·Rₖ·PTDFₖᵢ·Σⱼ PTDFₖⱼ·Pⱼ

where:
- Rₖ: resistance of arc k
- PTDFₖᵢ: power transfer distribution factor for arc k, bus i
- Pⱼ: injection at bus j

This is a simplified formulation that assumes constant voltage magnitudes
(typically 1.0 p.u.). It's computationally faster but less accurate than
AC loss factors for systems with significant voltage variations.
"""
function get_DC_dLoss_dP(res, sys)
    ptdf = PTDF(sys)
    bus_length = length(axes(ptdf, 1))
    arcs_length = length(axes(ptdf, 2))
    injection = get_bus_injection(res)
    T_length = size(injection, 2)
    loss_p = zeros(bus_length, T_length)
    
    R, _ = get_RX_vector(sys, ptdf)
    
    # Compute DC loss sensitivity for each bus and time period
    # Assumes constant voltage (V = 1.0 p.u.)
    for i in 1:bus_length
        for t in 1:T_length
            loss_p[i, t] = sum(2 * R[k] * ptdf[k, i] * (sum(ptdf[k, j] * injection[j, t] for j in 1:bus_length)) for k in 1:arcs_length)
        end
    end
    return loss_p
end

"""
    get_total_AC_loss(res) -> Vector{Float64}

Compute total system active power losses from AC power flow results.

# Arguments
- `res`: SimulationResults or DecisionModelResults with power flow variables

# Returns
- Vector of length num_timesteps containing total losses (in MW) at each time

# Details
Total losses are computed as the sum of losses on all branches:
- Lines: Loss = P_from_to + P_to_from (both directions)
- Tap transformers: Loss = P_from_to + P_to_from

For a lossless branch, P_from_to = -P_to_from, so the sum is zero.
For lossy branches, the sum is positive and represents energy dissipated.

This is the "ground truth" loss from detailed AC power flow calculations.
"""
function get_total_AC_loss(res)
    # Read power flows in both directions for lines
    FromTo_Line = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerFromTo__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_Line = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerToFrom__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    
    # Read power flows in both directions for tap transformers
    FromTo_TapTransformer = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerFromTo__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_TapTransformer = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerToFrom__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    
    T_length = size(FromTo_Line, 1)
    total_loss = zeros(T_length)
    
    # Sum all branch losses: Loss_k = P_from_to + P_to_from
    for t = 1:T_length
        total_loss[t] = sum(FromTo_Line[t, :] + ToFrom_Line[t, :]) + sum(FromTo_TapTransformer[t, :] + ToFrom_TapTransformer[t, :])
    end
    return total_loss
end

function get_total_AC_loss(res::PSI.SimulationProblemResults{PowerSimulations.DecisionModelSimulationResults})
    # Read power flows in both directions for lines
    FromTo_Line = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerFromTo__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_Line = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerToFrom__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    
    # Read power flows in both directions for tap transformers
    FromTo_TapTransformer = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerFromTo__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_TapTransformer = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerToFrom__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    
    T_length = size(FromTo_Line, 1)
    total_loss = zeros(T_length)
    
    # Sum all branch losses: Loss_k = P_from_to + P_to_from
    for t = 1:T_length
        total_loss[t] = sum(FromTo_Line[t, :] + ToFrom_Line[t, :]) + sum(FromTo_TapTransformer[t, :] + ToFrom_TapTransformer[t, :])
    end
    return total_loss
end

"""
    get_line_loss(res) -> Matrix{Float64}

Compute individual line losses from AC power flow results.

# Arguments
- `res`: SimulationResults or DecisionModelResults

# Returns
- Matrix of size (num_lines × num_timesteps) containing losses per line (in MW)
  - Row k: losses on line k
  - Column t: losses at time t

# Details
Loss on each line is computed as P_from_to + P_to_from.
For lossless lines, this sum is zero. For lossy lines, it represents
the I²R losses dissipated in the series resistance.
"""
function get_line_loss(res)    
    FromTo_Line = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerFromTo__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_Line = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerToFrom__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    # Transpose to get (lines × time) format
    return copy(transpose(FromTo_Line + ToFrom_Line))
end

function get_line_loss(res::PSI.SimulationProblemResults{PowerSimulations.DecisionModelSimulationResults})    
    FromTo_Line = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerFromTo__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_Line = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerToFrom__Line"; table_format = TableFormat.WIDE)[!, 2:end])
    # Transpose to get (lines × time) format
    return copy(transpose(FromTo_Line + ToFrom_Line))
end

"""
    get_tap_transformer_loss(res) -> Matrix{Float64}

Compute individual tap transformer losses from AC power flow results.

# Arguments
- `res`: SimulationResults or DecisionModelResults

# Returns
- Matrix of size (num_transformers × num_timesteps) containing losses per transformer (in MW)
  - Row k: losses on transformer k
  - Column t: losses at time t

# Details
Loss on each transformer is computed as P_from_to + P_to_from.
Transformers have both resistive (I²R) and core losses, though this
formulation primarily captures the load-dependent resistive losses.
"""
function get_tap_transformer_loss(res)    
    FromTo_TapTransformer = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerFromTo__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_TapTransformer = Matrix{Float64}(read_aux_variable(res, "PowerFlowLineActivePowerToFrom__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    # Transpose to get (transformers × time) format
    return copy(transpose(FromTo_TapTransformer + ToFrom_TapTransformer))
end

function get_tap_transformer_loss(res::PSI.SimulationProblemResults{PowerSimulations.DecisionModelSimulationResults})    
    FromTo_TapTransformer = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerFromTo__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    ToFrom_TapTransformer = Matrix{Float64}(read_realized_aux_variable(res, "PowerFlowLineActivePowerToFrom__TapTransformer"; table_format = TableFormat.WIDE)[!, 2:end])
    # Transpose to get (transformers × time) format
    return copy(transpose(FromTo_TapTransformer + ToFrom_TapTransformer))
end

"""
    get_bus_from_index(bus_numbers, line) -> Int

Find the index of the source bus in a bus numbering array.

# Arguments
- `bus_numbers`: Vector of bus numbers in the system
- `line`: ACBranch component (Line or Transformer)

# Returns
- Index position of the "from" bus in the bus_numbers array

# Details
Extract the source bus number from the branch's arc and locates its
position in the provided bus numbering scheme.
"""
function get_bus_from_index(bus_numbers, line)
    bus_from = line.arc.from
    bus_from_no = get_number(bus_from)
    return findfirst(x -> x == bus_from_no, bus_numbers)
end

"""
    get_bus_to_index(bus_numbers, line) -> Int

Find the index of the destination bus in a bus numbering array.

# Arguments
- `bus_numbers`: Vector of bus numbers in the system
- `line`: ACBranch component (Line or Transformer)

# Returns
- Index position of the "to" bus in the bus_numbers array

# Details
Extract the destination bus number from the branch's arc and locates its
position in the provided bus numbering scheme.
"""
function get_bus_to_index(bus_numbers, line)
    bus_to = line.arc.to
    bus_to_no = get_number(bus_to)
    return findfirst(x -> x == bus_to_no, bus_numbers)
end

"""
    get_fictitious_nodal_demand_by_loss(res, sys) -> Matrix{Float64}

Distribute transmission losses to buses as fictitious nodal demands.

# Arguments
- `res`: SimulationResults or DecisionModelResults
- `sys`: PowerSystems.jl System object

# Returns
- Matrix of size (num_buses × num_timesteps) containing fictitious demands (in MW)
  - FND[i,t] = fictitious demand at bus i, time t

# Details
This method allocates branch losses to the terminal buses for loss accounting:
- Each branch loss is split equally (50/50) between its two endpoint buses
- Includes losses from both transmission lines and tap transformers
- Only considers available (in-service) branches

Fictitious Nodal Demand (FND) is a bookkeeping mechanism where transmission
losses are represented as additional loads at buses. This approach allows
simpler formulations to approximate losses by adjusting nodal demands rather
than explicitly modeling losses in the power flow equations.

Use cases:
- Loss allocation for economic analysis
- Simplified loss representation in optimization models
- Post-processing of power flow results
"""
function get_fictitious_nodal_demand_by_loss(res, sys)
    # Read power flow variables
    line_var = read_aux_variable(res, "PowerFlowLineActivePowerFromTo__Line"; table_format = TableFormat.WIDE)
    tap_var = read_aux_variable(res, "PowerFlowLineActivePowerFromTo__TapTransformer"; table_format = TableFormat.WIDE)
    T_length = size(line_var, 1)
    
    # Extract component names
    line_names = names(line_var[!, 2:end])
    tap_names = names(tap_var[!, 2:end])
    
    # Get bus numbering
    bus_numbers = parse.(Int, names(read_expression(res, "ActivePowerBalance__ACBus", table_format=TableFormat.WIDE)[!, 2:end]))
    
    # Compute individual branch losses
    line_loss = get_line_loss(res)
    tap_tap_tx_loss = get_tap_transformer_loss(res)
    
    # Initialize fictitious nodal demand matrix
    FND = zeros(length(bus_numbers), T_length)
    
    # Allocate line losses (50% to each endpoint)
    for (ix_line, line_name) in enumerate(line_names)
        line = get_component(Line, sys, line_name)
        if !get_available(line)
            continue  # Skip unavailable lines
        end
        bus_ix_from = get_bus_from_index(bus_numbers, line)
        bus_ix_to = get_bus_to_index(bus_numbers, line)
        
        # Split loss equally between both buses
        for t = 1:T_length
            FND[bus_ix_from, t] += line_loss[ix_line, t] / 2.0
            FND[bus_ix_to, t] += line_loss[ix_line, t] / 2.0
        end
    end
    
    # Allocate transformer losses (50% to each endpoint)
    for (ix_tap, tap_name) in enumerate(tap_names)
        tap = get_component(TapTransformer, sys, tap_name)
        if !get_available(tap)
            continue  # Skip unavailable transformers
        end
        bus_ix_from = get_bus_from_index(bus_numbers, tap)
        bus_ix_to = get_bus_to_index(bus_numbers, tap)
        
        # Split loss equally between both buses
        for t = 1:T_length
            FND[bus_ix_from, t] += tap_tap_tx_loss[ix_tap, t] / 2.0
            FND[bus_ix_to, t] += tap_tap_tx_loss[ix_tap, t] / 2.0
        end
    end
    
    return FND
end

function get_fictitious_nodal_demand_by_loss(res::PSI.SimulationProblemResults{PowerSimulations.DecisionModelSimulationResults}, sys)
    # Read power flow variables
    line_var = read_realized_aux_variable(res, "PowerFlowLineActivePowerFromTo__Line"; table_format = TableFormat.WIDE)
    tap_var = read_realized_aux_variable(res, "PowerFlowLineActivePowerFromTo__TapTransformer"; table_format = TableFormat.WIDE)
    T_length = size(line_var, 1)
    
    # Extract component names
    line_names = names(line_var[!, 2:end])
    tap_names = names(tap_var[!, 2:end])
    
    # Get bus numbering
    bus_numbers = parse.(Int, names(read_realized_expression(res, "ActivePowerBalance__ACBus", table_format=TableFormat.WIDE)[!, 2:end]))
    
    # Compute individual branch losses
    line_loss = get_line_loss(res)
    tap_tap_tx_loss = get_tap_transformer_loss(res)
    
    # Initialize fictitious nodal demand matrix
    FND = zeros(length(bus_numbers), T_length)
    
    # Allocate line losses (50% to each endpoint)
    for (ix_line, line_name) in enumerate(line_names)
        line = get_component(Line, sys, line_name)
        if !get_available(line)
            continue  # Skip unavailable lines
        end
        bus_ix_from = get_bus_from_index(bus_numbers, line)
        bus_ix_to = get_bus_to_index(bus_numbers, line)
        
        # Split loss equally between both buses
        for t = 1:T_length
            FND[bus_ix_from, t] += line_loss[ix_line, t] / 2.0
            FND[bus_ix_to, t] += line_loss[ix_line, t] / 2.0
        end
    end
    
    # Allocate transformer losses (50% to each endpoint)
    for (ix_tap, tap_name) in enumerate(tap_names)
        tap = get_component(TapTransformer, sys, tap_name)
        if !get_available(tap)
            continue  # Skip unavailable transformers
        end
        bus_ix_from = get_bus_from_index(bus_numbers, tap)
        bus_ix_to = get_bus_to_index(bus_numbers, tap)
        
        # Split loss equally between both buses
        for t = 1:T_length
            FND[bus_ix_from, t] += tap_tap_tx_loss[ix_tap, t] / 2.0
            FND[bus_ix_to, t] += tap_tap_tx_loss[ix_tap, t] / 2.0
        end
    end
    
    return FND
end

"""
    remove_double_circuit_name(branch_name) -> String

Remove the "-double_circuit" suffix from a branch name.

# Arguments
- `branch_name`: String name of a branch component

# Returns
- Cleaned branch name with "-double_circuit" suffix removed

# Details
Some branch names include a "-double_circuit" suffix to indicate parallel circuits.
This function normalizes the names by removing this suffix, which is useful for
matching branches across different data representations where the suffix may or
may not be present.

# Examples
```julia
remove_double_circuit_name("Line123-double_circuit")  # Returns: "Line123"
remove_double_circuit_name("Transformer5")            # Returns: "Transformer5"
```
"""
function remove_double_circuit_name(branch_name)
    # Remove "-double_circuit" suffix if present
    return replace(branch_name, "-double_circuit" => "")
end

"""
    get_arc_axis_from_branch_name(sys, branch_name) -> Tuple{Int, Int}

Extract the (from_bus, to_bus) arc definition for a branch by name.

# Arguments
- `sys`: PowerSystems.jl System object
- `branch_name`: String name of the branch (may include "-double_circuit" suffix)

# Returns
- Tuple (from_bus_number, to_bus_number) representing the arc

# Details
This function:
1. Removes any "-double_circuit" suffix from the branch name
2. Searches for an ACBranch component whose name contains the trimmed name
3. Extracts and returns the bus numbers defining the branch's arc

Useful for mapping constraint keys (which use branch names) to PTDF matrix
indices (which use arc tuples). Handles naming inconsistencies between
different data sources or representations.

# Example
```julia
arc = get_arc_axis_from_branch_name(sys, "Line_A_B-double_circuit")
# Returns: (bus_A_number, bus_B_number)
```
"""
function get_arc_axis_from_branch_name(sys, branch_name)
    # Remove any double circuit suffix for matching
    trim_branch_name = remove_double_circuit_name(branch_name)
    
    # Find the branch component by partial name match
    branch = first(PSY.get_components(x -> contains(x.name, trim_branch_name), ACBranch, sys))
    
    # Return the arc as (from_bus_number, to_bus_number)
    return (branch.arc.from.number, branch.arc.to.number)
end

"""
    compute_iterative_error_based_on_generator_output(
        res_old,
        res_new
    ) -> Float64

Compute convergence metric for iterative loss approximation algorithms.

# Arguments
- `res_old`: Results from previous iteration
- `res_new`: Results from current iteration

# Returns
- Maximum infinity norm of changes across all generator types (in MW)

# Details
This function quantifies how much generator dispatch changed between iterations
by computing the infinity norm (maximum absolute change) for:
- Thermal generators (ThermalStandard)
- Renewable generators (RenewableDispatch)

The error metric is:
    error = max(‖P_thermal_new - P_thermal_old‖_∞, ‖P_renewable_new - P_renewable_old‖_∞)

where ‖·‖_∞ is the maximum absolute difference across all generators and time periods.

**Convergence Interpretation:**
- Small error (< 1e-3 MW): Dispatch has stabilized, loss approximation converged
- Large error: Significant redispatch occurred, need more iterations
- Increasing error: Possible divergence, check model formulation

**Why This Metric?**
Generator output changes indicate whether the loss approximation has stabilized.
When losses are accurately represented, the optimal dispatch should remain
consistent between iterations.

# Example
```julia
error = compute_iterative_error_based_on_generator_output(res_iter1, res_iter2)
if error < 1e-3
    println("Converged!")
end
```
"""
function compute_iterative_error_based_on_generator_output(res_old, res_new)
    # Read thermal generator outputs from previous iteration
    th_res = read_variable(res_old, "ActivePowerVariable__ThermalStandard"; table_format = TableFormat.WIDE)
    re_res = read_variable(res_old, "ActivePowerVariable__RenewableDispatch"; table_format = TableFormat.WIDE)
    th_old = Matrix{Float64}(th_res[!, 2:end])
    re_old = Matrix{Float64}(re_res[!, 2:end])
    
    # Read thermal generator outputs from current iteration
    th_res = read_variable(res_new, "ActivePowerVariable__ThermalStandard"; table_format = TableFormat.WIDE)
    re_res = read_variable(res_new, "ActivePowerVariable__RenewableDispatch"; table_format = TableFormat.WIDE)
    th_new = Matrix{Float64}(th_res[!, 2:end])
    re_new = Matrix{Float64}(re_res[!, 2:end])

    # Compute infinity norm (max absolute change) for each generator type
    err_th = norm(sum(th_new, dims = 2) - sum(th_old, dims = 2), Inf)  # Max change in thermal dispatch
    err_re = norm(sum(re_new, dims = 2) - sum(re_old, dims = 2), Inf)  # Max change in renewable dispatch
    
    # Return the maximum error across all generator types
    return maximum([err_th, err_re])
end