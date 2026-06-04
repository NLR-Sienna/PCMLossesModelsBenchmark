"""
    add_internal_hvdc!(sys::PSY.System)

Add two synthetic `TwoTerminalGenericHVDCLine` components to `sys` in-place.

The two links are Newark↔NRS (bus 8335↔8814) and Metcalf↔SanJoseB (bus 1819↔1258).
Both are configured with ±1000 MW (±10 p.u.) active-power limits and a lossless
`LinearCurve(0.0)` loss model.

These HVDC lines create a closed loop in the CATS network that can carry circular
flows when the optimization finds it cost-optimal to route power around the loop.

# Arguments
- `sys::PSY.System`: the system to modify in-place.

# Returns
Nothing — mutates `sys`.
"""
function add_internal_hvdc!(sys::PSY.System)
    newark_bus = PSY.get_bus(sys, 8335) # Newark Bus: Closest to [37.500801, -121.985703]
    nrs_bus = PSY.get_bus(sys, 8814) # NRS Bus: Closest to [37.402, -121.970001]
    newark_nrs_link = PSY.TwoTerminalGenericHVDCLine(
        name = "Newark_NRS_HVDC",
        available = true,
        active_power_flow = 0.0,
        arc = PSY.Arc(newark_bus, nrs_bus),
        active_power_limits_from = (-10.0, 10.0), # 1000 MW
        active_power_limits_to = (-10.0, 10.0), # 1000 MW
        reactive_power_limits_from = (-5.0, 5.0), # 500 MVar
        reactive_power_limits_to = (-5.0, 5.0), # 500 MVar
        loss = PSY.LinearCurve(0.0), # no loss        
    )
    PSY.add_component!(sys, newark_nrs_link)
    metcalf_bus = PSY.get_bus(sys, 1819) # Metcalf Bus: Closest to [37.224201, -121.741898]
    sanjoseb_bus = PSY.get_bus(sys, 1258) # SanJoseB Bus: Closest to [37.3409, -121.901604]
    metcalf_sanjoseb_link = PSY.TwoTerminalGenericHVDCLine(
        name = "Metcalf_SanJoseB_HVDC",
        available = true,
        active_power_flow = 0.0,
        arc = PSY.Arc(metcalf_bus, sanjoseb_bus),
        active_power_limits_from = (-10.0, 10.0), # 1000 MW
        active_power_limits_to = (-10.0, 10.0), # 1000 MW
        reactive_power_limits_from = (-5.0, 5.0), # 500 MVar
        reactive_power_limits_to = (-5.0, 5.0), # 500 MVar
        loss = PSY.LinearCurve(0.0), # no loss        
    )
    PSY.add_component!(sys, metcalf_sanjoseb_link)
    return
end