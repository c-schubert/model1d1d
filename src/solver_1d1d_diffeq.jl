using DifferentialEquations

# TODOS:
# macro generate dTdt_s! for MORE speed...
# clean up and doc ...

"""
EssentialModelComprehension:
Essential Properties of the different modell structs, needed to setup the solution
procedure
"""
struct EssentialModelComprehension{T1 <: AbstractFloat, T2 <: Integer}
    n_c::T2
    n_c_s_sub::T2
    T_in_f::T1
    T_in_s::T1
    dx_f::T1
    mf_f::T1
    mf_s::T1
    A_fs::T1
    A_f::T1
    htc_fs::HeatTransfer
    mp_f::IncompressibleFluidThermalMaterial
    mp_s::SolidThermalMaterial
    V_f::T1
    V_s::T1
    V_fs::T1
    scrap_cond::Bool
    dx_s_sub::T1
    A_s_sub::T1
    V_s_sub::T1
    cross_rad::Bool
    front_rad::Bool
    front_rad_model::Union{Nothing,HotSideRadiationModel}
    A_fs_rad::T1
    b::Vector{<:T1}
    k::T2
end


function solve_1d1d_diffeq(m::M1d1d)
    print_log(m1d1d_logfile, "Using diff eq solver!")

    # heat_flow_exchng_sf = zeros(Float64, m.mesh_settings.no_cells_bulk)
    emc = Vector{EssentialModelComprehension{Float64, Int64}}(undef, length(m.sim_phase))
    
    # Init Results according to used meshes
    results = Vector{Results1d1d}(undef,0)

    mesh_tmp = nothing

    phase_start_time = 0.0
    for pid=1:1:length(m.sim_phase)
        print_log(m1d1d_logfile, ["\nStarting Simulation of Phase ", pid])
        local T

        mesh = M1d1dMesh(m, pid)
        printInfo(mesh; tabs = 1)

        # just make sure that meshes do not change between the phases
        # because we than would have to interpolate between the meshes if they
        # change:
        if pid > 1
            if mesh_tmp != mesh
                error("Mesh change between different phases is not supported atm")
                print_log(m1d1d_logfile, "Error: Unsupported mesh change")
            else
                mesh_tmp = mesh
            end
        else
            mesh_tmp = mesh
        end

        emc[pid] = set_emc(m, pid, mesh)

        if pid > 1 
            phase_start_time += m.sim_phase[pid-1].duration
            T = set_T(emc[pid], pid, results[end])
        else
            T = set_T(emc[pid], pid)
        end
        
        qalpha_ex = zeros(emc[pid].n_c)
        dTdt = zeros(size(T))
        
        tspan = (phase_start_time, phase_start_time + m.sim_phase[pid].duration)

        # p = [emc[pid], qalpha_ex]
        range = zeros(Int, 2*emc[pid].k)
        p = (emc[pid], qalpha_ex, emc[pid].mp_f.lambda, emc[pid].mp_f.cp_true, 
                emc[pid].mp_f.rho, emc[pid].mp_s.lambda, emc[pid].mp_s.cp_true, 
                emc[pid].mp_s.rho, emc[pid].mp_s.epsilon, range)
        @show T[1,:]
        @show tspan
        prob = ODEProblem(dTdt_model!,T , tspan, p)
    
        @time sol = solve(prob,Tsit5(), progress = true, saveat=m.sim_phase[pid].dt_save
                , reltol=1e-8, abstol=1e-8,dtmin=1e-6, maxiters=1e8)

        # test perf / allocs ...
        # fix many allocations of material functions pX_bounden and integrate_pX
        # @time dTdt_model!(dTdt,T,p,0.2)
        # @time dTdt_model2!(dTdt,T,p,0.2)
        push_sol_to_res!(results, m, mesh, emc[pid].mp_s.cp_true, emc[pid].mp_s.rho, 
                emc[pid].V_s_sub, sol)
    end
  
    return results
 end


function push_sol_to_res!(res::Vector{Results1d1d},m::M1d1d, mesh::M1d1dMesh, 
        cp_s::MaterialPropertyOfX, rho_s::MaterialPropertyOfX, V_s_sub_c::Real, 
        sol)

    for t = 1:length(sol.t)
        Tres = M1d1dVals(mesh, init_val = 0.0)

        T_H = collect(0:0.05:2000)
        H_T_s = H_of_T.(cp_s, T_H)

        @inbounds for i = 1:1:size(sol.u[t],1)
            Tres.f[i] =  sol.u[t][i,1]

            if size(sol.u[t],2) > 2
                for j = 1:1:size(sol.u[t],2) 
                    Tres.s_sub[i][j] = sol.u[t][i,j+1]
                end

                # determine Tmean of all T_s_subs if we use only T_s at a later 
                # stage
                m_sum = 0.0
                H_T_sub_sum = 0.0
                @simd for j = 2:1:size(sol.u[t],2) 
                    m_sum += pX_bounded(rho_s,sol.u[t][i,j])  * V_s_sub_c
                    H_T_sub_sum += H_of_T(cp_s, sol.u[t][i,j]) * 
                            pX_bounded(rho_s, sol.u[t][i,j]) * V_s_sub_c
                end

                Tres.s[i] = T_of_H(H_T_s, T_H, H_T_sub_sum/(m_sum))
            else
                Tres.s[i] =  sol.u[t][i,2]

                if !isnothing(Tres.s_sub)
                # make sure that we can start with T_s_sub in later pid run 
                    for j = 1:1:length(Tres.s_sub[1]) 
                        Tres.s_sub[i][j] = sol.u[t][i,2]
                    end
                end
            end
        end
      
        push!(res, Results1d1d(sol.t[t], Tres))
    end
end


function set_T(emc::EssentialModelComprehension{<:AbstractFloat, <:Integer}, 
        phaseid::Integer, Tprev::Union{Nothing, Results1d1d} = nothing)

    if emc.scrap_cond 
        T = zeros(emc.n_c, emc.n_c_s_sub+1)
    else
        T = zeros(emc.n_c, 2) 
    end

    if phaseid == 1
        T .= emc.T_in_s

    elseif !isnothing(Tprev) 
        for  i = 1:emc.n_c
            T[i,1] = Tprev.T.f[i]
        end

        if emc.scrap_cond 
            for i = 1:emc.n_c, j=2:1:emc.n_c_s_sub+1
                T[i,j] = Tprev.T.s_sub[i][j-1]
            end
        else
            for  i = 1:emc.n_c
                T[i,2] = Tprev.T.s[i]
            end
        end
        
    else
        error("No Tprev given!")
    end

    return T
 end


function dTdt_model!(dTdt::Matrix{<:Real}, T::Matrix{<:Real}, p::Tuple{EssentialModelComprehension{<:AbstractFloat, <:Integer}, 
    Vector{<:AbstractFloat}, MaterialPropertyOfX, MaterialPropertyOfX,MaterialPropertyOfX, 
    MaterialPropertyOfX,MaterialPropertyOfX,MaterialPropertyOfX,MaterialPropertyOfX,Vector{<:Integer}}, t)
    emc = p[1]
    qalpha_ex = p[2]
    lam_f = p[3]
    cp_f = p[4]
    rho_f = p[5]
    lam_s = p[6]
    cp_s = p[7]
    rho_s = p[8]
    eps_s  = p[9]
    range = p[10]

    dTdt_f!(dTdt, T, qalpha_ex, emc, lam_f, cp_f, rho_f)
    dTdt_s!(dTdt, T, qalpha_ex, emc, lam_s, cp_s, rho_s, eps_s, range)
end


function dTdt_f!(dTdt::Matrix{<:AbstractFloat}, T::Matrix{<:AbstractFloat}, qalpha_ex::Vector{<:AbstractFloat}, emc::EssentialModelComprehension{<:AbstractFloat, <:Integer},
        lam_f::MaterialPropertyOfX, cp_f::MaterialPropertyOfX, rho_f::MaterialPropertyOfX)

    lambda_dT(TI, T2)::Float64 = harm_mean(pX_bounded(lam_f,T2), 
            pX_bounded(lam_f, TI)) * (T2 - TI)
    f_dHdt_f_mvmnt(T1, T2)::Float64 = ( emc.mf_f * integrate_pX(cp_f, T1, T2) )

    @inbounds @simd for i = 1:1:emc.n_c
        local dHdT_alpha::Float64
        local dHdT_cond::Float64
        local dHdt_f_mvmnt::Float64
        local dHdt_f::Float64
        local dTdt_f::Float64

        dHdT_alpha = 0.0
        dHdT_alpha = (emc.A_fs * htc_fs(emc.htc_fs) * 
                (T[i,2] - T[i,1]))

        qalpha_ex[i] = dHdT_alpha
        dHdT_cond = 0.0
        dHdt_f_mvmnt = 0.0

        if i != 1 && i != emc.n_c
            dHdT_cond = emc.A_f/emc.dx_f * ( lambda_dT(T[i,1], T[i-1,1]) +
                lambda_dT(T[i,1], T[i+1,1]) )

            dHdt_f_mvmnt = f_dHdt_f_mvmnt(T[i,1], T[i-1,1])
        elseif i == emc.n_c
            dHdT_cond = emc.A_f/emc.dx_f * lambda_dT(T[i,1], T[i-1,1]) 

            dHdt_f_mvmnt = f_dHdt_f_mvmnt(T[i,1], T[i-1,1])
        elseif i == 1
            dHdT_cond = emc.A_f/emc.dx_f * lambda_dT(T[i,1], T[i+1,1])
            
            dHdt_f_mvmnt = f_dHdt_f_mvmnt(T[i,1], emc.T_in_f) 
            pX_bounded(lam_f, 100.0) 
        end

        dHdt_f = (dHdt_f_mvmnt + dHdT_alpha + dHdT_cond)
        dTdt_f = 0.0
        dTdt_f = dHdt_f/(emc.V_f * pX_bounded(rho_f, T[i,1]) * 
                    pX_bounded(cp_f, T[i,1]))

        dTdt[i,1] = dTdt_f
    end
end


"""
we could optimize that further if we use a macro to write a specific version
for this function ..
"""
function dTdt_s!(dTdt::Matrix{T1}, T::Matrix{T1}, qalpha_ex::Vector{T1}, 
    emc::EssentialModelComprehension{T1, T2}, lam_s::MaterialPropertyOfX, 
    cp_s::MaterialPropertyOfX, rho_s::MaterialPropertyOfX, 
    eps_s::MaterialPropertyOfX, range::Vector{T2}) where T1 <: AbstractFloat where T2 <: Integer

    lambda_dT(TI, T2)::T1 = ( harm_mean(pX_bounded(lam_s,T2), 
    pX_bounded(lam_s, TI)) * (T2 - TI) )

    #  f(i,k) = vcat(collect(i-k:1:i-1), collect(i+1:1:i+k))
    function update_range!(range, i, k)
        @inbounds @simd for ii = 0:1:k
            range[ii] = (ii+i-1) - k
        end

        @inbounds @simd for ii = k+1:1:2*k
            range[ii] = (ii-k) + i  
        end
    end

    @inbounds for i = emc.n_c:-1:1
        dHdT_alpha = -qalpha_ex[i]

        local dHdt_rad_front::T1 
        local dHdt_rad_cross::T1
        local dHdt_cond::T1
        local dHdt_s_mvmnt::T1
        local dHdt_s::T1

        dHdt_rad_front = 0.0
        if emc.front_rad
            if i <= emc.k 
                if i == 1
                    dHdt_rad_front= heat_flow_front_rm(
                            emc.front_rad_model, 
                            pX_bounded(eps_s, T[i,2]), emc.b[i], 
                            emc.A_fs_rad, T[i,2]) 
                else
                    dHdt_rad_front = heat_flow_front_rm(
                            emc.front_rad_model, 
                            pX_bounded(eps_s, T[i,2]), emc.b[i+1], 
                            emc.A_fs_rad, T[i,2]) 
                end
            end
        end

        # cross radiation model
        dHdt_rad_cross = 0.0
        if emc.cross_rad
            update_range!(range, i, emc.k)

            for ir = 1:1:2*emc.k
                local r
                r = range[ir]
                if r > 0 && r <= emc.n_c
                    T_r = T[r,2]
                    T_i = T[i,2]

                    r_dist = abs(i-r)
                    
                    dHdt_rad_cross += emc.A_fs_rad *
                            hf_radiation_par_surf_b_weighted( 
                            pX_bounded(eps_s, T_r), pX_bounded(eps_s, T_i), 
                            emc.b[r_dist+1], T_r, T_i )
                end
            end
        end

        if emc.scrap_cond
            for j = 2:1:emc.n_c_s_sub+1
                dHdt_cond = 0.0
                scrap_cond_fac = emc.A_s_sub/emc.dx_s_sub
                if j != 2 && j != (emc.n_c_s_sub+1)
                    # CDS
                    dHdt_cond = scrap_cond_fac * ( lambda_dT(T[i,j], T[i,j-1]) +
                                lambda_dT(T[i,j], T[i,j+1]) )
                elseif j == (emc.n_c_s_sub+1) && j > 2
                    # symmetry only 1 dir (BDS + 0)
                        dHdt_cond = scrap_cond_fac * lambda_dT(T[i,j], T[i,j-1]) 
                elseif j == 2
                    # convection (FDS)
                    dHdt_cond = scrap_cond_fac * lambda_dT(T[i,j], T[i,j+1]) 
                end
        
                dHdt_s_mvmnt = 0.0
                if i != emc.n_c
                    dHdt_s_mvmnt = (emc.mf_s * integrate_pX(cp_s, 
                            T[i,j], T[i+1,j]))
                else
                    dHdt_s_mvmnt = (emc.mf_s * integrate_pX(cp_s, 
                            T[i,j], emc.T_in_s))
                end

                dHdt_s = 0.0
                if j == 2
                    dHdt_s = ( 0.5 * (dHdT_alpha + 
                            dHdt_rad_front + dHdt_rad_cross) + 
                            dHdt_cond + dHdt_s_mvmnt)
                else
                    dHdt_s = (dHdt_s_mvmnt + dHdt_cond)
                end
                
                dTdt_s = dHdt_s/(emc.V_s_sub * pX_bounded(rho_s, T[i,j]) * 
                        pX_bounded(cp_s, T[i,j]))

                dTdt[i,j] = dTdt_s
            end
        else
            dHdt_s_mvmnt = 0.0
            if i != emc.n_c
                dHdt_s_mvmnt = (emc.mf_s * integrate_pX(cp_s, 
                        T[i,2], T[i+1,2]))
            else
                dHdt_s_mvmnt = (emc.mf_s * integrate_pX(cp_s, 
                        T[i,2], emc.T_in_s))
            end

            dHdt_s = (dHdT_alpha +  dHdt_rad_front + dHdt_rad_cross  + dHdt_s_mvmnt)
            
            dTdt_s = dHdt_s/(emc.V_s * pX_bounded(rho_s, T[i,2]) * 
                    pX_bounded(cp_s, T[i,2]))

            dTdt[i,2] = dTdt_s
        end
    end
end 


function set_emc(m::M1d1d, phaseid::Integer, mesh::M1d1dMesh)
 
    scrap_conduction = false
    scrap_conduction_temp = false
    A_fs_rad = 0.0

    no_cells = mesh.no_cells
    dx_f = mesh.gas_cells.len

    A_f = mesh.gas_cells.dA
    A_s = mesh.solid_cells.dA

    V_f = mesh.gas_cells.vol
    V_s = mesh.solid_cells.vol
    V_fs = V_f + V_s

    mf_f = m.offgas[phaseid].mass_flow
    mf_s = m.bulk[phaseid].mass_flow

    max_b_fac_layers = floor(Int64, no_cells/2)
    b_fac_over_layers = zeros(max_b_fac_layers)
    k = 0

    print_log(m1d1d_logfile, ["Mass flow fluid: ", mf_f])
    print_log(m1d1d_logfile, ["Mass flow solid: ", mf_s])

    A_fs =  mesh.A_fs
    htc_fs = m.heat_transfer[phaseid]

    print_log(m1d1d_logfile, ["Af_s: ", A_fs]);

    no_sub_cells = 0
    dx_s_sub = 0.0
    A_s_sub =  0.0
    V_s_sub = 0.0 

    if m.sim_phase[phaseid].submodels.scrap_discretization_model_enabled
        print_log(m1d1d_logfile, "Scrap conduction enabled")

        if isnothing(mesh.solid_1d_submesh)
            error("No submesh generated for scrap part heating!")
        else
            no_sub_cells = mesh.solid_1d_submesh.no_cells
            dx_s_sub = mesh.solid_1d_submesh.sub_cells.len
            A_s_sub = mesh.solid_1d_submesh.sub_cells.dA
            V_s_sub = mesh.solid_1d_submesh.sub_cells.vol
            mf_s = mf_s / (no_sub_cells * 2.0)

            print_log(m1d1d_logfile, "Scrap conduction settings: ")
            print_log(m1d1d_logfile, ["\t no_sub_cells: ", no_sub_cells])
            print_log(m1d1d_logfile, ["\t dx_s_sub: ", dx_s_sub])
            print_log(m1d1d_logfile, ["\t A_s_sub: ", A_s_sub])
            print_log(m1d1d_logfile, ["\t A_s_sub/dx_s_sub: ", A_s_sub/dx_s_sub])
            print_log(m1d1d_logfile, ["\t V_s_sub: ", V_s_sub])
            print_log(m1d1d_logfile, ["\t mf_s: ", mf_s])

            if no_sub_cells < 3
                error("Not able to handle submesh with less than 3 cells")
            end
        end
    end
    
    print_log(m1d1d_logfile, "Radiation Settings: ")

    cross_radiation = m.sim_phase[phaseid].submodels.cross_radiation_model_enabled
 
    if cross_radiation
        print_log(m1d1d_logfile,"\t Cross radiation enabled")
    else
        print_log(m1d1d_logfile,"\t Cross radiation disabled")
    end

    front_radiation = !isnothing(
            m.sim_phase[phaseid].submodels.furnace_side_rad_model)

    if front_radiation
        print_log(m1d1d_logfile, "\t Front radiation enabled")
    else
        print_log(m1d1d_logfile, "\t Front radiation disabled")
    end

    if cross_radiation || front_radiation

        c = m.sim_phase[phaseid].submodels.cross_radiation_constant_c

        @assert c <= (1.0 + eps(typeof(c)))  && c >= (0.0 - eps(typeof(c)))

        A_fs_rad = c * 0.5 * A_fs

        A_cross = m.bulk[phaseid].cross_sect_surf_area
        eps_s = m.bulk[phaseid].scrap_props.thermal_props.epsilon

        if A_fs_rad > A_cross/pX_bounded(eps_s, Tmax)
            # TODO: we asume here pX_bounded(eps_s, Tmax) == epsilon_max
            # me may want to enhance the way this limit here is set ...
            warnmsg = ("A_fs_rad bigger than  A_cross/epsilon_max, this would lead " *
                    "to unphysical results, consider choosing a smaller " *
                    "cross radiation constant c. Limiting A_fs_rad to A_cross/epsilon_max")

            @warn   warnmsg

            print_log(m1d1d_logfile, "Warning: " * warnmsg, print_to_stdout=false)
            A_fs_rad =  A_cross/pX_bounded(eps_s, Tmax)

            # TODO: is that right because b[1] stands for self viewing? or not? ...
            k = 1
            b_fac_over_layers[1] = 1.0 
            b_fac_over_layers[2] = 1.0 

            b_fac_over_layers_print = b_fac_over_layers[b_fac_over_layers .> 1E-3]
            print_log(m1d1d_logfile, ["\t b_fac_over_layers: ", b_fac_over_layers_print])
            print_log(m1d1d_logfile, ["\t k: ",  k])

        elseif A_fs_rad > 0
            mean_overlap = calc_mean_overlap(A_cross, A_fs_rad, 1, 10_000) 

            for j = 1:length(b_fac_over_layers)
                if j == 1
                    b_fac_over_layers[j] = 1 # with itself
                elseif j == 2
                    b_fac_over_layers[j] = mean_overlap  # with direct neighbor
                else
                    b_fac_over_layers[j] = (1-mean_overlap).^(j-2) *
                            mean_overlap # thrugh cells bfore
                end
            end
            
            b_fac_over_layers_print = b_fac_over_layers[b_fac_over_layers .> 1E-3]
            print_log(m1d1d_logfile, ["\t b_fac_over_layers: ", b_fac_over_layers_print])
            
            k = findfirst( x -> x < 1E-3, b_fac_over_layers)
            print_log(m1d1d_logfile, ["\t k: ",  k])
        else
            error("Wrong settings of cross_radiation_constant_c?")
        end

        if isnothing(k)

            warnmsg = ("b_fac_over_layers[" * string(max_b_fac_layers) * 
                    "] is very high, maybe there is " * 
                    "something wrong? k will be set to maximum of " *
                    string(max_b_fac_layers))

            @warn warnmsg
            print_log(m1d1d_logfile, "Warning: " * warnmsg, print_to_stdout=false)

            k = max_b_fac_layers
        end

        print_log(m1d1d_logfile, ["\t Af_s_rad: ", A_fs_rad]);
    end

    return EssentialModelComprehension(no_cells, no_sub_cells, 
            m.offgas[phaseid].T_in, m.bulk[phaseid].T_charged, dx_f, mf_f, mf_s, 
            A_f, A_fs, htc_fs, m.offgas[phaseid].thermal_props, 
            m.bulk[phaseid].scrap_props.thermal_props, V_f, V_s, V_fs, 
            m.sim_phase[phaseid].submodels.scrap_discretization_model_enabled,
            dx_s_sub, A_s_sub, V_s_sub,
            cross_radiation, front_radiation, 
            m.sim_phase[phaseid].submodels.furnace_side_rad_model,
            A_fs_rad, b_fac_over_layers, k)
end
