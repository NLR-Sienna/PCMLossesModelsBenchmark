using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder
using PowerFlows
using PowerNetworkMatrices
using HydroPowerSimulations
using InfrastructureSystems
using JuMP
using LinearAlgebra
using HiGHS
using Dates
using Logging
import PowerSystems as PSY
import PowerSimulations as PSI
import PowerSystemCaseBuilder as PSB

# the rating is assumed to be 100 MVA for all lines, and the line length is 25 km for all lines. The voltage level is assumed to be 345 kV for all lines.
# use miso transmission cost estimate AC new single circuit transmission line cost (table 4.1.1) in Texas
function build_expansion_ac_line_model(sys, line_arc, line_length, voltage_level, line_name)
    if voltage_level == 230
        base_voltage = 230
        # use the 703MVA line
        line_r = 5.63327E-05 * line_length
        line_x = 0.000638941 * line_length
        line_b = 0.002575331 * line_length
        line_rating = 10
        cost = 2.6 *1000000 * line_length /(8760*40) #2.6
    elseif voltage_level == 345
        base_voltage = 345
        # use the 1010MWA line
        line_r = 2.68011E-05 * line_length
        line_x = 0.000314304 * line_length
        line_b = 0.00526269 * line_length
        line_rating = 10
        cost = 4 * 1000000 * line_length /(8760*40)
    elseif voltage_level == 500
        base_voltage = 500
        # use the 2503MWA line
        line_r = 0.00000832 * line_length
        line_x = 0.00011776 * line_length
        line_b = 0.014004 * line_length
        line_rating = 10
        cost = 5.1 * 1000000 * line_length /(8760*40) #5.1
    elseif voltage_level == 765
        base_voltage = 765
        #  use the 5300MVA line
        line_r = 2.18719E-06 * line_length
        line_x = 4.85967E-05 * line_length
        line_b = 0.034394844 * line_length
        line_rating = 10
        cost = 6.3 *1000000 * line_length /(8760*40)
    else
        error("Unsupported voltage level: $voltage_level kV")
    end
    line_name = line_name * "_" * string(voltage_level) * "kV"
    candidate_line = Line(
            name = line_name,
            available = true,
            active_power_flow = 0,
            reactive_power_flow = 0,
            arc = line_arc,
            r = line_r,
            x = line_x,
            b = (from = line_b, to = line_b),
            rating = line_rating,
            angle_limits = (min = -0.5235987755982988, max = 0.5235987755982988),
            ext = Dict(
                "is_candidate" => true,
                "project_cost" => cost,
            ),
        ) 
        return candidate_line
end

