# include("../activate_env.jl")

using JLD2
using model1d1d

#= ###########################################################################
Example for Simulation Setup
########################################################################### =#
casename = "my_example_case"

sdm_enabled = false
front_rad_enabled = true 
cross_rad_enabled = true

logpath="./"
logfile = joinpath(logpath, casename*"_log.txt")
new_log(logfile) 

desc = casename * ", logfile: ", logfile


submodel = model1d1d.SimSettingsSubModels(
        scrap_discretization_model_enabled = false,
        cross_radiation_model_enabled = true, 
        furnace_side_rad_model =  model1d1d.SimpleRadiationModel(900+273.15)
        )

gas_scrap_htc = model1d1d.ConstantHeatTransfer(50.0) # in W/(m² K)

#= ###########################################################################
Process Settings
########################################################################### =#
shaft_length  = 7.0 # in m 
shaft_width   = 3.0 # in m 
shaft_height  = 3.0 # in m 

shaft_cross_sect_area = shaft_width * shaft_height

duration = 50.0 * 60.0 # in seconds (48 min) - old 40 / 2

# save results every t_save seconds
t_save = 30.0

# setting general simulation settings
sim_phases = [model1d1d.SimSettings(duration, submodel, t_save)]
        
#= ###########################################################################
Hot-Gas Settings
########################################################################### =#
gas_volume_flow = 19500 # NCm_per_h

gas_inlet_temperature =  1200 + 273.15 # K

# Material data
gas_density = model1d1d.ReciprocalMaterialPropertyOfX(360.47962061763553, 0.0,
        name = "density", unit = "kg/m^3")

gas_cp_true = model1d1d.PolynomialMaterialPropertyOfX([908.692796686028,
    0.336368996611928, -6.950621118017843e-5], name = "true heat capacity", 
    unit = "J/(kg K)", bounds = (1000.0,1500.0))

gas_lambda = model1d1d.PolynomialMaterialPropertyOfX([0.007187645849107943, 
        5.784085171400682e-5, 5.428742980927024e-9, -4.200871592280629e-12], 
        name = "heat conductivity", unit = "W/(m K)", bounds = (0.01,0.15))

gas_mat = model1d1d.IncompressibleFluidThermalMaterial("off-gas", gas_density,
        gas_cp_true, gas_lambda)

# mass flow for norm temperature
gas_mass_flow = gas_volume_flow/3600 * model1d1d.pX(gas_mat.rho, 25 + 273.15)

offgas = model1d1d.OffGas(gas_mass_flow, gas_inlet_temperature, gas_mat)

#= ###########################################################################
Bulk Settings
########################################################################### =#
# complete mass of charged scrap
scrap_mass_charged = 111.1 * 1000.0 # in kg

scrap_bulk_mass_flow = scrap_mass_charged / duration

scrap_charged_temperature = 25.0 + 273.15 # in K

scrap_bulk_density = 801.0	# kg / m³ at RT
scrap_avg_thickness_of_piece = 0.01 # in m
scrap_avg_length_of_piece = 0.5 # in m
scrap_avg_width_of_piece = 0.4 # in m

# Material data
scrap_density = model1d1d.ConstantMaterialPropertyOfX(7850.0, name = "density", 
        unit ="kg/m^3")

scrap_cp_true = model1d1d.ConstantMaterialPropertyOfX(500.0,  name = "true heat capacity",
unit = "J/(kg K)")

scrap_lambda = model1d1d.ConstantMaterialPropertyOfX(30.0, name = "heat conductivity"
, unit = "W/(m K)")

scrap_epsilon = model1d1d.ConstantMaterialPropertyOfX(0.8, name = "emissivity", 
        unit ="kg/m^3")

scrap_mat = model1d1d.SolidThermalMaterial("scrap", scrap_density, 
        scrap_cp_true, scrap_lambda, scrap_epsilon)


print_log(logfile, "\nBulk properties:")
print_log(logfile, ["\t Scrap bulk density: ", scrap_bulk_density])
print_log(logfile, ["\t Scrap avgerage thickness of piece: ", scrap_avg_thickness_of_piece])
print_log(logfile, ["\t Scrap avgerage length of piece: ", scrap_avg_length_of_piece])
print_log(logfile, ["\t Scrap avgerage width of piece: ", scrap_avg_width_of_piece])

scrap_bulk_porosity = 1.0 - (scrap_bulk_density/model1d1d.pX(scrap_mat.rho, 25 + 273.15))

print_log(logfile, ["\t Scrap bulk porosity: ", scrap_bulk_porosity])

scrap_piece = model1d1d.CharScrapPiece(scrap_avg_thickness_of_piece,
        scrap_avg_length_of_piece, scrap_avg_width_of_piece, scrap_mat)

bulk = model1d1d.Bulk(scrap_bulk_mass_flow, 
        scrap_charged_temperature, shaft_length, shaft_width*shaft_height,
        scrap_bulk_density, scrap_bulk_porosity, 0.95, scrap_piece)

#= ###########################################################################
Numerical Settings
########################################################################### =#

no_cells_over_shaft = 120
no_cells_over_scrap_piece = nothing


mesh_settings = model1d1d.M1d1dMeshSettings(no_cells_over_shaft, no_cells_over_scrap_piece)

#= ###########################################################################
Final Model
########################################################################### =#

model = model1d1d.M1d1d(sim_phases, [gas_scrap_htc], [bulk], 
        [offgas], [mesh_settings]);


#= ###########################################################################
Solution
########################################################################### =#
results = model1d1d.solve_1d1d_diffeq(model);


#= ###########################################################################
Save Results
########################################################################### =#

@save joinpath("./results", casename*".jld2") model desc results 


