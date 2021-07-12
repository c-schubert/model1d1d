#= #############################################################################
Mesh Types

The mesh is generated via preprocessing from the M1d1d struct just before the 
simulation.
############################################################################# =#
abstract type AbstractM1d1dMesh end;

"""
M1dCell:
Represents the structure for a regular, equaly sized 1d mesh:

    vol: volume of cell
    len: length of cell 
    dA: cross section area normal to direction of len
"""
struct M1dCell{T <: AbstractFloat}
    vol::T
    len::T 
    dA::T
end

"""
SubMesh1D: 
Mesh type for the additiona discretization of the scrap pieces
here a symmetry conditin will be used, so the real no of cells 
is ceil(no_cells_scrap_piece/2)
"""
struct SubMesh1D
    no_cells::Integer
    sub_cells::M1dCell
end

"""
Bulk1dMesh:
Mesh consiting out of gas_cells and solid_cells which overlap each other and may
exchange energy terms during the modeling. Both cells types define 1d 
discretized meshes with constant spacing and cell volume.

    no_cells: number cells over length of shaft (bulk)
    A_fs: Interface area between solid and fluid cell(s)
    gas_cells: Properties of fluid (gas) cells
    solid_cells: Properties of solid cells
    solid_1d_submesh: If existent (!= nothing) proerties of the solid 1d submesh for thermal heat transfer inside the scrap pieces in each solid bulk (shaft) cell
"""
struct M1d1dMesh <: AbstractM1d1dMesh
    no_cells::Integer
    A_fs::AbstractFloat
    gas_cells::M1dCell
    solid_cells::M1dCell
    solid_1d_submesh::Union{Nothing,SubMesh1D}
end


#= #############################################################################
Mesh Value Setting Types
############################################################################# =#
"""
M1d1dVals:
Structure for storing (property) values for each M1d1d cell
    f: values of fluid cells
    s: values of solid cells
    s_sub: values for solid sub cells
"""
struct M1d1dVals{T <: AbstractFloat}
    f::Vector{T}
    s::Vector{T}
    s_sub::Union{Nothing, Vector{Vector{T}}}
end

Base.copy(S::M1d1dVals) = M1d1dVals(S.f, S.s, S.s_sub)


#= #############################################################################
Results Types
############################################################################# =#

"""
ResultsTimeStepWithoutScrapDisc:
Results for temperature distribution in bulk and offgas
"""
struct Results1d1d 
    t::Float64
    T::M1d1dVals
end


#= #############################################################################
Outer Constructor Functions
############################################################################# =#
function M1d1dMesh(m::M1d1d, pid::Integer)

    # number bulk cells
    no_c = m.mesh_settings[pid].no_cells_bulk

    # length of bulk cell
    len_c =  m.bulk[pid].length / no_c

    # bulk volume
    V_b = m.bulk[pid].length * m.bulk[pid].cross_sect_surf_area

    # bulk volume of single cell
    V_bc = V_b / no_c

    V_sc = (1.0 - m.bulk[pid].porosity) * V_bc

    V_fc = m.bulk[pid].porosity * V_bc

    dA_sc = (1-m.bulk[pid].porosity) *  m.bulk[pid].cross_sect_surf_area

    dA_fc = m.bulk[pid].porosity * m.bulk[pid].cross_sect_surf_area

    solid_cells = model1d1d.M1dCell(V_sc, len_c, dA_sc)

    gas_cells = model1d1d.M1dCell(V_fc, len_c, dA_fc)

    weight_scrap_piece =  (m.bulk[pid].scrap_props.vol * pX_bounded(
        m.bulk[pid].scrap_props.thermal_props.rho, m.bulk[pid].T_charged))

    weight_bulk_cell = V_bc * m.bulk[pid].density 

    no_scrap_pieces_per_bulk_cell = weight_bulk_cell / weight_scrap_piece

    # surface area between solid and gas cells
    A_fs_c =  ( m.bulk[pid].scrap_pieces_contact_factor * 
            no_scrap_pieces_per_bulk_cell * m.bulk[pid].scrap_props.surf_area )

    if !isnothing(m.mesh_settings[pid].no_cells_scrap_piece)
        # cross section surface in thickness direction for solid scrap cell(s)
        no_sub_c = m.mesh_settings[pid].no_cells_scrap_piece

        dA_s_sub = m.bulk[pid].scrap_props.char_cross_section_area * 
                no_scrap_pieces_per_bulk_cell

        len_s_sub = m.bulk[pid].scrap_props.char_thickness / (2.0 * no_sub_c)

        vol_s_sub = V_sc / (2.0 * no_sub_c)

        solid_1d_subcell = M1dCell(vol_s_sub, len_s_sub, dA_s_sub)

        solid_1d_submesh = SubMesh1D(no_sub_c, solid_1d_subcell)

    else
        solid_1d_submesh = nothing
    end
    
    return M1d1dMesh(no_c, A_fs_c, gas_cells, solid_cells, solid_1d_submesh)
end


function M1d1dVals(mesh1d1d::M1d1dMesh;
    init_val = 0.0)::M1d1dVals

    init_vec = ones(Float64, mesh1d1d.no_cells) .* init_val

    if !isnothing(mesh1d1d.solid_1d_submesh)
        sub_mesh_init = Vector{Vector{Float64}}(undef, mesh1d1d.no_cells)

        for i = 1:length(sub_mesh_init)
            sub_mesh_init[i] = (ones(Float64, mesh1d1d.solid_1d_submesh.no_cells) 
            .* init_val)
        end
    else
        sub_mesh_init = nothing  
    end

    return M1d1dVals(copy(init_vec), copy(init_vec), sub_mesh_init)
end

#= #############################################################################
Other functions
############################################################################# =#

function Base.isequal(m1::M1d1dMesh, m2::M1d1dMesh)

    return
    (
        m1.no_cells ==  m2.no_cells &&
        isapprox(m1.A_fs, m2.A_fs) &&
        m1.gas_cells == m2.gas_cells &&
        m1.soid_cells == m2.soid_cells
    )

end

function Base.isequal(sm1::Union{Nothing, SubMesh1D}, sm2::Union{Nothing, SubMesh1D})
    return
    (
        (isnothing(sm1) && isnothing(sm2)) ||
        (
            sm1.no_cells == sm2.no_cells &&
            sm1.sub_cells == sm2.sub_cells
        )
    )
end


function Base.isequal(c1::M1dCell, c2::M1dCell)
    return
    (
        isapprox(c1.vol, c2.vol) &&
        isapprox(c1.len, c2.len) &&
        isapprox(c1.dA, c2.dA)
    )
end


"""
print1d1dvals:
Formated output of M1d1dVals struct.
"""
function print1d1dvals(v::M1d1dVals)

    println("T-off-gas: \t T-solid")
    for i = 1:1:length(v.f)
        Printf.@printf("% 6.2f \t % 6.2f",v.f[i],v.s[i])

        if !isnothing(v.s_sub)
            print("---->(" )
            for j=1:1:length(v.s_sub[i])
                Printf.@printf("% 6.2f,",v.s_sub[i][j])
            end
            print(")" )
        end
        print("\n")
    end
end


function printInfo(mesh::M1d1dMesh;tabs::Int=0)

    if tabs > 0
        tabs_str = "\t".^tabs
    else
        tabs_str = ""
    end

    print_log(m1d1d_logfile, tabs_str * "Mesh-Info:")

    print_log(m1d1d_logfile, tabs_str * "No bulk cells: " * string(mesh.no_cells))
    print_log(m1d1d_logfile, tabs_str * "Cell gas-solid surface area: " * string(mesh.A_fs))
    print_log(m1d1d_logfile, tabs_str * "Gas cells:")
    printInfo(mesh.gas_cells;tabs = (tabs + 1))
    print_log(m1d1d_logfile, tabs_str * "Solid cells:")
    printInfo(mesh.solid_cells;tabs = (tabs + 1))

    if !isnothing(mesh.solid_1d_submesh)
        print_log(m1d1d_logfile, tabs_str * "Submesh:")
        print_log(m1d1d_logfile, tabs_str * "\t no sub cells: " *
                string(mesh.solid_1d_submesh.no_cells))
            print_log(m1d1d_logfile, tabs_str * "\t Subcells:")
            printInfo(mesh.solid_1d_submesh.sub_cells;tabs = (tabs + 2))
    end
end


function printInfo(mesh_cell::M1dCell; tabs::Int=0)

    if tabs > 0
        tabs_str = "\t".^tabs
    else
        tabs_str = ""
    end
    print_log(m1d1d_logfile, tabs_str * "vol: " * string(mesh_cell.vol) )
    print_log(m1d1d_logfile, tabs_str * "len: " * string(mesh_cell.len) )
    print_log(m1d1d_logfile, tabs_str * "dA: " * string(mesh_cell.dA) )
end

