"""
HotSideRadiationModel:
Implementation of types to model the incident radiation on the "hot" side of 
the scrap bulk.
"""
abstract type HotSideRadiationModel end


"""
SimpleRadiationModel - For a simple radiation model:

dot{Q}_rad = sigma * epsilon * A * (T - T_inf)

"""
struct SimpleRadiationModel <: HotSideRadiationModel
    T_inf::Float64
end

### S2S Model types

"""
S2Ssurface:
Represents a general surface in the s2s model.
"""
struct S2Ssurface{T1 <: Int, T2 <: Real}
    name::String
    id::T1
    T::T2
    area::T2
    epsilon::T2
end

# """
# S2SRadiationModel
# To use a more sophisticated radiation model, where externally calculated view 
# factors can be used, to model the influence an envoirement with different hot 
# spots in front of the scrap bulk.
# """
mutable struct S2SRadiationModel{T1 <: Integer, T2 <: Real} <: HotSideRadiationModel
    front_surface_id::T1
    n_surfaces::T1
    surfaces::Vector{S2Ssurface{T1, T2}}
    vf_mat::Matrix{T2}
    has_initialized::Bool
    A::Matrix{T2}
    B::Vector{T2}
    # G::Vector{T2}
    # Qrad::Vector{T2}
end


function set_S2SRadiationModel(surf_names::Vector{String},  vf_mat::Matrix{T1},
T_surfaces::Vector{T1}, A_surfaces::Vector{T1}, epsilon_surfaces::Vector{T1},
front_surface_id::T2 = 1) where {T1 <: Real, T2 <: Integer}

    @assert size(vf_mat,1) == size(vf_mat,2) 

    n = size(vf_mat,1)

    @assert n == length(surf_names)
    @assert n == length(T_surfaces)
    @assert n == length(A_surfaces)
    @assert n == length(epsilon_surfaces)

    surfaces = Vector{S2Ssurface{T2, T1}}(undef, n)
    for i = 1:n
        surfaces[i] = S2Ssurface(surf_names[i], i, T_surfaces[i], A_surfaces[i],
                epsilon_surfaces[i])
    end

    A = zeros(n, n)
    B = zeros(n)

    return S2SRadiationModel(front_surface_id, n, surfaces, vf_mat, false, A, B) # , G, Qrad)
end


"""
heat_flow_rad_front(rm::S2SRadiationModel,eps::Real, A::Real,T::Real)

heat flow calculation for scrap front cells with s2s modelß0
"""
function heat_flow_front_rm(rm::S2SRadiationModel, eps::Real, b::Real, A::Real, T::Real)::Real

    if !rm.has_initialized
        @inbounds for i1 = 1:rm.n_surfaces, i2 = 1:rm.n_surfaces
            if i1 == i2
                rm.A[i1, i2] = ((rm.surfaces[i1].epsilon - 1) * rm.vf_mat[i1, i2]) + 1
            else
                rm.A[i1, i2] = (rm.surfaces[i1].epsilon - 1) * rm.vf_mat[i1, i2]
            end
        end

        T_faces = zeros(rm.n_surfaces)

        for i = 1:rm.n_surfaces
            T_faces[i]  = rm.surfaces[i].T
        end

        T_faces[rm.front_surface_id] = T

        for ii = 1:rm.n_surfaces
            rm.B[ii] = sigma * rm.surfaces[ii].epsilon * (T_faces[ii]^4);
        end

        rm.has_initialized = true
    else
        rm.B[rm.front_surface_id] = sigma * rm.surfaces[rm.front_surface_id].epsilon * (T^4);
    end

    J = rm.A\rm.B;

    # what do we realy need for Q_other_to_front_surface_id
    # @inbounds for i1 = 1:rm.n_surfaces
    #     G_zw = 0
    #     for i2 = 1:rm.n_surfaces
    #         G_zw = G_zw + J[i2] * rm.vf_mat[i1, i2]
    #     end
    #     rm.G[i1] = G_zw
    # end

    # @inbounds for ii=1:rm.n_surfaces
    #     rm.Qrad[ii] = (J[ii] - rm.G[ii]) * rm.surfaces[ii].area
    # end

    # Qrad = -1 .* rm.Qrad[rm.front_surface_id] 
    
    # G: incident radiation
    # J: outgoining including reflection

    G_front_face = 0.0
    for i2 = 1:rm.n_surfaces
        G_front_face = G_front_face + J[i2] * rm.vf_mat[rm.front_surface_id, i2]
    end

    Qrad = (J[rm.front_surface_id] - G_front_face) * rm.surfaces[rm.front_surface_id].area

    return - b * Qrad
end
                                  

"""
heat_flow_rad_front(rm::SimpleRadiationModel,eps::Real, A::Real,T::Real)

heat flow calculation for scrap front cells ...
"""
function heat_flow_front_rm(rm::SimpleRadiationModel, eps::Real, b::Real, A::Real, T::Real)::Real

    return b * hf_radiation_inf(eps, A, rm.T_inf, T)
end

