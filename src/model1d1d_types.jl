"""
Types used to build the structure for the 1d1d scrap bulk modeling 
"""


"""
    M1d1dMeshSettings

Settings for the mesh of the 1d1d modelling. Mesh is generated in the simulation step, from bulk, scrap and shaft properties itself.

# Arguments
- `no_cells_bulk::Integer`: number of cells over bulk length
- `no_cells_scrap_piece::Union{Nothing, Integer}`: number of cells over scrap piece
"""
struct M1d1dMeshSettings
    no_cells_bulk::Integer
    no_cells_scrap_piece::Union{Nothing, Integer}
end

#= #############################################################################
Simulation Setting Types
############################################################################# =#

"""
    SimSettingsSubModels
    
Settings for activation and special options for the submodels

"""
struct SimSettingsSubModels{T <: AbstractFloat}
    scrap_discretization_model_enabled::Bool    # defaults to false
    min_no_of_scrap_cells::Integer              # defaults to 3
    cross_radiation_model_enabled::Bool         # defaults to false
    cross_radiation_constant_c::T               # defaults to ?
    furnace_side_rad_model::Union{Nothing,HotSideRadiationModel}

    SimSettingsSubModels(;scrap_discretization_model_enabled = false,
            min_no_of_scrap_cells = 5,
            cross_radiation_model_enabled = false, 
            cross_radiation_constant_c = 1.0, 
            furnace_side_rad_model = nothing) = 
            new{typeof(cross_radiation_constant_c)}(scrap_discretization_model_enabled,
                    min_no_of_scrap_cells, cross_radiation_model_enabled, 
                    cross_radiation_constant_c, furnace_side_rad_model)
end

"""
SimSettings:
Struct representing the simulation settings for a phase in the 
simulation of the 1d1d model.

All uinits musst be given in SI units, temperatures in K

dt_save:
interval time for results beeing written to Vector{ResultsTimeStep**ScrapDisc}
"""
struct SimSettings{T <: AbstractFloat}
    duration::T
    submodels::SimSettingsSubModels{T}
    dt_save::T # defaults  to duration?

    SimSettings(
                duration::T2, 
                submodels::SimSettingsSubModels{T2}
                ) where {T2 <: AbstractFloat}  = new{T2}(duration, submodels, duration)

    SimSettings(
        duration::T2, submodels::SimSettingsSubModels{T2}, dt_save::T2
        ) where {T2 <: AbstractFloat} = new{T2}(duration, submodels, dt_save)
end



#= #############################################################################
Scrap Bulk Types
############################################################################# =#

# """
# ScrapMixture:
# Desribes the mixture of different ScrapPiece(s). Here mass_frac[i] is the mass 
# fraction of a ScrapPiece (pieces[i]).
# (sum(mass_frac[:] musst be equal to 1)
#TODO: Make this work in the future
# """
# mutable struct ScrapMixture
#     pieces::Vector{ScrapPiece}
#     mass_frac::Vector{Float64}
# end


"""
CharScrapPiece:
A mutable struct CharScrapPiece (charcterisitic piece of scrap) is asumed 
to be a flat shaped object which has a charcterisitic 
thickness that can represent the heating throughout it's body.
Also it has a surf_area which ist mostly relevant for convection and 
radiation transport.
"""
struct CharScrapPiece{T <: AbstractFloat}
    char_thickness::T
    char_cross_section_area::T
    thermal_props::SolidThermalMaterial
    surf_area::T
    vol::T

    # initialize from dimensions and mat
    CharScrapPiece( T::T2,
                    W::T2,
                    L::T2,
                    mat::SolidThermalMaterial) where {T2<:AbstractFloat}=
                new{T2}(
                        char_thickness_quader([T, W, L]),
                        char_cross_section_area_quader([T, W, L]),
                        mat,
                        surface_quader(T,W,L),
                        volume_quader(T,W,L)
                    )
end

"""
Bulk:
Struct representing the properties and geometry of the scrap bulk

scrap_pieces_contact_factor: factor which defines the average amount of scrap 
piece surface, which has direct contact with another scrap piece.
All units musst be given in SI units
"""
struct Bulk{T <: AbstractFloat}
    mass_flow::T
    T_charged::T
    length::T
    cross_sect_surf_area::T
    density::T
    porosity::T
    scrap_pieces_contact_factor::T
   # mat::ScrapMixture # Future see ScrapMixture
    scrap_props::CharScrapPiece
end

#= #############################################################################
Off-gas Types
############################################################################# =#

"""
OffGas:
Struct representing the off-gas

All uinits musst be given in SI units, temperatures in K
"""
struct OffGas{T <: AbstractFloat}
    mass_flow::T
    T_in::T
    thermal_props::IncompressibleFluidThermalMaterial
end


#= #############################################################################
Final Model Type
############################################################################# =#

"""
M1d1d:
Model parent type
"""
struct M1d1d
    sim_phase::Vector{SimSettings}
    heat_transfer::Vector{HeatTransfer}
    bulk::Vector{Bulk}
    offgas::Vector{OffGas}
    mesh_settings::Vector{M1d1dMeshSettings}
end
