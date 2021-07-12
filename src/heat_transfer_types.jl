#= #############################################################################
Heat Transfer Types
############################################################################# =#

"""
HeatTransfer:
Types for the implementation of different heat transfer correlations

TODO:
Currently we neglect influence of the material properties, since we would need 
information on the nusselt number (Nu = f(Re,Pr)).So maybe a nusselt number 
correlation could be determined with variable coefficients could be derived to 
model the heat transfer in scrap bulk.
"""
abstract type HeatTransfer end

"""
HeatTransferCurveFromVal:
The HeatTransfer coefficient through the scrap usally dependend thermal 
properties and the flow speed of the gas (turbulence). We only consider the 
flow speed here, by presribing the alpha value (val) accoding to a certain flow 
velocity. 
Also as for now we define a single parameter (mod_fac) to lineary shift this 
curve for manual optimization.
"""
struct HeatTransferCurveFromVal <: HeatTransfer
    val::Vector{Float64}
    flow_speed::Vector{Float64}
    mod_fac::Float64
end

"""
ConstantHeatTransfer:
Constant heat transfer coefficient
"""
struct ConstantHeatTransfer <: HeatTransfer
    val::Float64
end
