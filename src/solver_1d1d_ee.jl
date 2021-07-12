"""
Functions to solve the 1d1d scrap - off-gas model
Here we use an explicit euler (EE), which is not very elegant or fast, but allows
easy modifications atm - and should deliver accurate (physical) solution if stable ...

DISCLAIMER: 
The newer solver_1d1d_diffeq solver should be used, which uses DifferentialEquations.jl 
with the Tsit5() solver and should be generally faster and more robust.
This solver is only left here for development and comparision.

Be careful when using it, there may be instabillities and some tweaking of timestep 
sizes may be necessary to get correct results ...
"""


function solve1d1d_ee(m::M1d1d)

    print_log(m1d1d_logfile, "Using explicit euler solver!")
    # Init Results according to used meshes
    results = Vector{Results1d1d}(undef,0)

    mesh_tmp = nothing
    
    T_H = collect(0:0.05:Tmax)

    scrap_conduction = false
    scrap_conduction_temp = false
    phase_start_time = 0.0

    mesh = M1d1dMesh(m, 1)
    T = M1d1dVals(mesh, init_val = m.bulk[1].T_charged)
    Tnew = M1d1dVals(mesh, init_val = m.bulk[1].T_charged)  

    for p=1:1:length(m.sim_phase)
        t_save = m.sim_phase[p].dt_save

        mesh = M1d1dMesh(m, p)
        printInfo(mesh; tabs = 1)
    
        # just make sure that meshes do not change between the phases
        # because we than would have to interpolate between the meshes if they
        # change:
        if p > 1
            if mesh_tmp != mesh
                error("Mesh change between different phases is not supported atm")
                print_log(m1d1d_logfile, "Error: Unsupported mesh change")
            else
                mesh_tmp = mesh
            end
        else
            mesh_tmp = mesh
        end

        no_cells = mesh.no_cells
        dx_f = mesh.solid_cells.len
    
        A_f = mesh.gas_cells.dA
        A_s = mesh.solid_cells.dA
    
        V_f = mesh.gas_cells.vol
        V_s = mesh.solid_cells.vol
        V_fs = V_f + V_s
    
        if !isnothing(mesh.solid_1d_submesh)
            no_sub_cells = mesh.solid_1d_submesh.no_cells
            dx_s = mesh.solid_1d_submesh.sub_cells.len
            A_s_s =  mesh.solid_1d_submesh.sub_cells.dA
            V_s_sub_c = mesh.solid_1d_submesh.sub_cells.vol        
        end
    
        heat_flow_exchng_sf = zeros(Float64, m.mesh_settings.no_cells_bulk)

        println("Starting Simulation of Phase $p")

        mf_f = m.offgas[p].mass_flow
        mf_s = m.bulk[p].mass_flow
        # TODO Check solid and fluid are equal

        println("Mass flow fluid: ", mf_f)
        println("Mass flow solid: ", mf_s)

        @assert isapprox(V_fs, dx * m.bulk[p].cross_sect_surf_area)
        @assert  m.bulk[p].scrap_pieces_contact_factor <= 1.0 + 
                eps(typeof(m.bulk[p].scrap_pieces_contact_factor))

        # Fluid to Solid Surface (and vice versa)
        no_scrap_pieces_per_vol::Float64 = ( m.bulk[p].density / 
                (m.bulk[p].scrap_props.vol * pX_bounded(
                m.bulk[p].scrap_props.thermal_props.rho, m.bulk[p].T_charged)) )
        surface_of_scrap_per_vol::Float64 = (
                m.bulk[p].scrap_pieces_contact_factor * no_scrap_pieces_per_vol 
                * m.bulk[p].scrap_props.surf_area )

        A_fs =  surface_of_scrap_per_vol * V_fs

        lam_f =  m.offgas[p].thermal_props.lambda
        rho_f = m.offgas[p].thermal_props.rho
        cp_f = m.offgas[p].thermal_props.cp_true

        lam_s = m.bulk[p].scrap_props.thermal_props.lambda
        rho_s = m.bulk[p].scrap_props.thermal_props.rho
        cp_s = m.bulk[p].scrap_props.thermal_props.cp_true
        
        # working but actually slows down computation ...
        # println("Tabulation fluid material properties")
        # lam_f = TabulatedMaterialPropertyOfX(m.offgas[p].thermal_props.lambda, Trange)
        # rho_f = TabulatedMaterialPropertyOfX(m.offgas[p].thermal_props.rho, Trange)
        # cp_f = TabulatedMaterialPropertyOfX(m.offgas[p].thermal_props.cp_true, Trange)
       
        # lam_s = TabulatedMaterialPropertyOfX(m.bulk[p].scrap_props.thermal_props.lambda, Trange)
        # rho_s = TabulatedMaterialPropertyOfX(m.bulk[p].scrap_props.thermal_props.rho, Trange)
        # cp_s = TabulatedMaterialPropertyOfX(m.bulk[p].scrap_props.thermal_props.cp_true, Trange)
       
        alpha_fs = m.heat_transfer[p]

        println("Integrating H of T curve")
        H_T_f = H_of_T.(cp_f, T_H)
        H_T_s = H_of_T.(cp_s, T_H)

        # just to be sure ...
        check_monolithic_growing(H_T_f)
        check_monolithic_growing(H_T_s)

        println("Check Radiation Settings ...")

        if m.sim_phase[p].submodels.scrap_discretization_model_enabled
            println("Scrap conduction enabled")

            if isnothing(mesh.solid_1d_submesh)
                error("No submesh generated for scrap part heating!")
            end

            scrap_conduction = true
            # 2.0 due to symmetry
            mf_sub_s::Float64  = mf_s / (mesh.solid_1d_submesh.no_cells * 2.0)
        else
            scrap_conduction = false
        end

        if p > 1 && scrap_conduction && !scrap_conduction_temp 
            println("Interpolating scrap temperatures to sub cells")
            # set prev temperature of T.s[i] to all T.s_sub[i][j] 
            @inbounds for i = 1:no_cells
                for j = 1:no_sub_cells
                    T.s_sub[i][j] = T.s[i]
                end
            end
        end

        if p > 1 
            scrap_conduction_temp = scrap_conduction
            phase_start_time += m.sim_phase[p-1].duration
        end

        cross_radiation = m.sim_phase[p].submodels.cross_radiation_model_enabled

        if cross_radiation
            println("Cross radiation enabled")
        end

        front_radiation = !isnothing(
                m.sim_phase[p].submodels.furnace_side_rad_model)

        if front_radiation
            println("Front radiation enabled")
        end

        if cross_radiation || front_radiation

            c = m.sim_phase[p].submodels.cross_radiation_constant_c

            @assert c <= (1.0 + eps(typeof(c)))  && c >= (0.0 - eps(typeof(c)))

            A_fs_rad = c * 0.5 * A_fs

            A_cross = m.bulk[p].cross_sect_surf_area

            eps_s = m.bulk[p].scrap_props.thermal_props.epsilon

            max_b_fac_layers = floor(Int64, no_cells/2)
            b_fac_over_layers = zeros(max_b_fac_layers)

            if A_fs_rad > A_cross/pX_bounded(eps_s, Tmax)
                # TODO: we asume here pX_bounded(eps_s, Tmax) == epsilon_max
                # me may want to enhance the way this limit here is set ...

                @warn   "A_fs_rad bigger than  A_cross/epsilon_max, this would lead "
                        "to unphysical results, consider choosing a smaller " *
                        "cross radiation constant c. Limiting A_fs_rad to A_cross/epsilon_max"
                A_fs_rad =  A_cross/pX_bounded(eps_s, Tmax)
                k = 1
                
                b_fac_over_layers[1] = 1
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
                
                println("b_fac_over_layers: ", b_fac_over_layers[b_fac_over_layers .> 1E-3])
                
                k = findfirst( x -> x < 1E-3, b_fac_over_layers)
                println("k: ",  k)
            else
                error("Wrong settings of cross_radiation_constant_c?")
            end

            if isnothing(k)
                @warn "b_fac_over_layers[" * string(max_b_fac_layers) * 
                        "] is very high, maybe there is " * 
                        "something wrong? k will be set to maximum of " *
                        string(max_b_fac_layers)
                k = max_b_fac_layers
            end
        end
        
        progress = nothing 
        if !isapprox(m.sim_phase[p].duration, 0)
            progress = Progress(ceil(Int64, m.sim_phase[p].duration), 1) 
        end
        t::Float64 = 0.0
        
        dt = set_timestepsize(T, m , p, mesh)
        n_time_iters = 0
        println("Initial timestep size: ", dt)
        ran_seconds = 0

        while t <= m.sim_phase[p].duration

            if !isnothing(progress) && ceil(Int64,t) > ran_seconds
                update!(progress, ceil(Int64,t))
                ran_seconds = ceil(Int64,t)
            end
            
            dt = set_timestepsize(T, m , p, mesh)

            @inbounds for i = 1:1:no_cells
                # heat exchange fluid - solid (gas - scrap)
                heat_flow_exchng_fs = 0.0
                if !scrap_conduction
                    heat_flow_exchng_fs = (A_fs * htc_fs(alpha_fs) * 
                            (T.s[i] - T.f[i]))
                else
                    heat_flow_exchng_fs = (A_fs * htc_fs(alpha_fs) * 
                            (T.s_sub[i][1] - T.f[i]))
                end

                lambda_dT(TI, T2)::Float64 = harm_mean(pX_bounded(lam_f,T2), 
                        pX_bounded(lam_f, TI)) * (T2 - TI)

                # heat conduction through fluid
                heat_flow_cond_f = 0.0
                if i != 1 && i != no_cells
                    heat_flow_cond_f = A_f/dx_f * ( lambda_dT(T.f[i], T.f[i-1]) +
                        lambda_dT(T.f[i], T.f[i+1]) )
                elseif i == no_cells
                    heat_flow_cond_f = A_f/dx_f * lambda_dT(T.f[i], T.f[i-1]) 
                elseif i == 1
                    heat_flow_cond_f = A_f/dx_f * lambda_dT(T.f[i], T.f[i+1]) 
                end

                # movement of fluid (gas)
                #integration only for polynomial cp atm, see matprop_functons.jl
                f_dHdt_f_mvmnt(T2)::Float64 = ( mf_f * 
                        integrate_pX(cp_f, T.f[i], T2) )
                
                dHdt_f_mvmnt = 0.0

                if i != 1
                    dHdt_f_mvmnt = f_dHdt_f_mvmnt(T.f[i-1])
                else
                    dHdt_f_mvmnt = f_dHdt_f_mvmnt(m.offgas[p].T_in)  
                end
                
                dH_f::Float64 = dt * (dHdt_f_mvmnt + heat_flow_exchng_fs + 
                        heat_flow_cond_f)

                # get the new T from the changed enthalpy of the cell
                Tnew.f[i] = T_of_H( H_T_f, T_H, ( H_of_T(cp_f, T.f[i]) + 
                        dH_f/(V_f * pX_bounded(rho_f, T.f[i])) ) )

                heat_flow_exchng_sf[i] = - heat_flow_exchng_fs
                
                if Tnew.s[i] > Twarn
                    @warn ("Temperatures getting to high possible numerical error")
                end 
            end

            @inbounds for i = no_cells:-1:1
                    # radiation models:
                    local heat_flow_rad_front::Float64 
                    local heat_flow_rad_in_cross::Float64
                    # front radiation model 
                    heat_flow_rad_front = 0.0
                    if front_radiation
                        if i <= k 
                            if scrap_conduction
                                T_front = T.s_sub[i][1]
                            else
                                T_front = T.s[i]
                            end

                            if i == 1
                                heat_flow_rad_front= heat_flow_front_rm(
                                        m.sim_phase[p].submodels.furnace_side_rad_model, 
                                        pX_bounded(eps_s, T_front), b_fac_over_layers[i], 
                                        A_fs_rad, T_front) 
                            else
                                heat_flow_rad_front = heat_flow_front_rm(
                                        m.sim_phase[p].submodels.furnace_side_rad_model, 
                                        pX_bounded(eps_s, T_front), b_fac_over_layers[i+1], 
                                        A_fs_rad, T_front) 
                            end
                        end
                    end

                    # cross radiation model
                    heat_flow_rad_in_cross = 0.0
                    if cross_radiation
                        range1 = vcat(collect(i-k:1:i-1), collect(i+1:1:i+k))
          
                        for r in range1
                            if r > 0 && r <= no_cells
                                if scrap_conduction
                                    T_r = T.s_sub[r][1]
                                    T_i = T.s_sub[i][1]
                                else
                                    T_r = T.s[r]
                                    T_i = T.s[i]
                                end

                                r_dist = abs(i-r)
                                
                                heat_flow_rad_in_cross += 
                                        hf_radiation_par_surf_b_weighted( 
                                        pX_bounded(eps_s, T_r), pX_bounded(eps_s, T_i), 
                                        b_fac_over_layers[r_dist+1], T_r, T_i )
                            end
                        end
                    end

                if scrap_conduction
                    lambda_dT(TI, T2)::Float64 = ( harm_mean(pX_bounded(lam_s,T2), 
                            pX_bounded(lam_s, TI)) * (T2 - TI) )

                    # scrap conduction
                    for j = 1:1:no_sub_cells
                        heat_flow_sub_cond_s = 0.0
                        if j != 1 && j != no_sub_cells
                            heat_flow_sub_cond_s = A_s_s/dx_s * 
                                    ( lambda_dT(T.s_sub[i][j], T.s_sub[i][j-1]) +
                                      lambda_dT(T.s_sub[i][j], T.s_sub[i][j+1]) )
                        elseif j == no_sub_cells
                            # symmetry
                                heat_flow_sub_cond_s = A_s_s/dx_s * 
                                        lambda_dT(T.s_sub[i][j], T.s_sub[i][j-1]) 
                        elseif j == 1
                            # convection
                            heat_flow_sub_cond_s = A_s_s/dx_s * 
                                    lambda_dT(T.s_sub[i][j], 
                                    T.s_sub[i][j+1]) 
                        end
                        
                        local dHdt_s_mvmnt::Float64
                        dHdt_s_mvmnt = 0.0
                        #scrap movement
                        if i != no_cells
                            dHdt_s_mvmnt = (mf_sub_s * integrate_pX(cp_s, 
                                    T.s_sub[i][j], T.s_sub[i+1][j]))
                        else
                            dHdt_s_mvmnt = (mf_sub_s * integrate_pX(cp_s, 
                                    T.s_sub[i][j], m.bulk[p].T_charged))
                        end
                        
                        local dH_s::Float64
                        dH_s = 0.0
                        if j == 1
                            # 0.5 because symmetry
                            dH_s = dt * ( 0.5 * (heat_flow_exchng_sf[i] + 
                                    heat_flow_rad_front + heat_flow_rad_in_cross) + 
                                    heat_flow_sub_cond_s + dHdt_s_mvmnt)
                        else
                            dH_s = dt * (dHdt_s_mvmnt + heat_flow_sub_cond_s)
                        end

                        Tnew.s_sub[i][j] = T_of_H(H_T_s, T_H,
                                ( H_of_T(cp_s, T.s_sub[i][j]) + 
                                dH_s/(V_s_sub_c * pX_bounded(rho_s, T.s_sub[i][j])) ) )

                        if Tnew.s_sub[i][j] > Twarn
                            @warn ("Temperatures getting to high possible numerical error")
                        end 
                    end
                else
                    #integration only for polynomial cp atm, see matprop_functons.jl
                    f_dHdt_s_mvmnt(T2)::Float64 = (mf_s * integrate_pX(cp_s, T.s[i], T2))

                    local dHdt_s_mvmnt::Float64
                    dHdt_s_mvmnt = 0.0
                    #scrap movement
                    if i != mesh.no_cells
                        dHdt_s_mvmnt = f_dHdt_s_mvmnt(T.s[i+1])
                    else
                        dHdt_s_mvmnt = f_dHdt_s_mvmnt(m.bulk[p].T_charged)
                    end

                    # final new temperature
                    dH_s = dt * ( dHdt_s_mvmnt + heat_flow_exchng_sf[i] 
                            + heat_flow_rad_front + heat_flow_rad_in_cross )

                    Tnew.s[i] = T_of_H(H_T_s, T_H, (H_of_T(cp_s, T.s[i]) + 
                            dH_s/(V_s * pX_bounded(rho_s, T.s[i]))) )

                    if Tnew.s[i] > Twarn
                        @warn ("Temperatures getting to high possible numerical error")
                    end 
                end
            end

            # T.f[:] = copy(Tnew.f)
            unsafe_update_v1_to_v2!(T.f, Tnew.f, no_cells)

            if scrap_conduction
                # T.s_sub[:] = deepcopy(Tnew.s_sub)
                unsafe_update_vv1_to_vv2!(T.s_sub, Tnew.s_sub, no_cells, no_sub_cells )
                
                # get average scrap cell temperature for post prossing ...
                @inbounds for i = 1:no_cells
                    # usually factor 2 because of symmetry but its cuts out 
                    m_sum = 0.0
                    H_T_sub_sum = 0.0
                    @simd for j=1:no_sub_cells
                        m_sum += pX_bounded(rho_s,T.s_sub[i][j])  * V_s_sub_c
                        H_T_sub_sum += H_of_T(cp_s, T.s_sub[i][j]) * 
                                pX_bounded(rho_s, T.s_sub[i][j]) * V_s_sub_c
                    end

                    T.s[i] = T_of_H(H_T_s, T_H, H_T_sub_sum/(m_sum))
                end
            else
                # T.s[:] = copy(Tnew.s)
                unsafe_update_v1_to_v2!(T.s, Tnew.s, no_cells)
            end

            t += dt
            
            if t >= t_save || t > m.sim_phase[p].duration
            push!(results, Results1d1d(phase_start_time + t, deepcopy(T)))
                t_save = t + m.sim_phase[p].dt_save
            end

            n_time_iters += 1
        end

        println("Finsihed p ", p, " after ", n_time_iters, " timesteps!")
    end

    return results
end


"""
set_timestepsize:
set (hoefully) stable according to CFL and Fo Criteria
TODO: Radiation may lead to instabilities, we could to iterative time step Criteria
but  that takes also it time and as EE is not the best after all, we should implement
better method in the future, by now hardcode stable timestep factor ...
"""
@inline function dt_cfl_min(rho::MaterialPropertyOfX, c_vol::Real, m_flow::Real, temp_vec::Vector{<:Real})

    dt_min = 10_000_000.0 # big number
    if !isapprox(m_flow, 0)
        for temp in temp_vec
            local dt 
            dt = c_vol * pX_bounded(rho, temp)/m_flow
            if (dt_min > dt)
                dt_min = dt
            end
        end
    end 

    return dt_min
end

@inline function dt_fo_min(rho::T1, lam::T2, cp::T3, dx::Real, temp_vec::Vector{<:Real}) where  {T1,T2,T3 <: MaterialPropertyOfX}
    dt_min = 10_000_000.0 # big number
    for temp in temp_vec
        local dt 
        dt = (0.5 * pX_bounded(rho, temp) * pX_bounded(cp, temp) * dx^2)/
                pX_bounded(lam, temp)
        if (dt_min > dt)
            dt_min = dt
        end
    end
    return dt_min
end

@inline function set_timestepsize(T::M1d1dVals, m::M1d1d, p::Int64, mesh::M1d1dMesh)::Float64

    lam_f =  m.offgas[p].thermal_props.lambda
    rho_f = m.offgas[p].thermal_props.rho
    cp_f = m.offgas[p].thermal_props.cp_true

    lam_s = m.bulk[p].scrap_props.thermal_props.lambda
    rho_s = m.bulk[p].scrap_props.thermal_props.rho
    cp_s = m.bulk[p].scrap_props.thermal_props.cp_true

    # CFL & FO like timestep size ...
    dtmin_1 = dt_cfl_min(rho_f, mesh.gas_cells.vol, m.offgas[p].mass_flow, T.f)
    dtmin_2 = dt_cfl_min(rho_s, mesh.solid_cells.vol, m.bulk[p].mass_flow, T.s)

    dtmin_3 = dt_fo_min(rho_f, lam_f, cp_f, mesh.gas_cells.len, T.s)

    dtmin_4 = 10_000_000.0 # big number
    if m.sim_phase[p].submodels.scrap_discretization_model_enabled
        T_s_min = minimum(T.s_sub[end])
        T_s_max = maximum(T.s_sub[1])

        if T_s_min < T_s_max
            steps = 5.0
            if (T_s_max - T_s_min)/steps > 1.0
                T_s = collect(T_s_min:steps:T_s_max)
            else
                T_s = [(T_s_max + T_s_min)/2.0]
            end 
        else
            # equal?
            T_s = [T_s_max]
        end

        dtmin_4 = dt_fo_min(rho_s, lam_s, cp_s, mesh.solid_1d_submesh.sub_cells.len, T_s)
    end

    dt = minimum((dtmin_1, dtmin_2, dtmin_3, dtmin_4))
    
    if !isnothing(m.sim_phase[p].submodels.furnace_side_rad_model)
        # hardcoded due to EE instabilities with radiation ...
        dt = 0.05 * dt
        if m.sim_phase[p].submodels.scrap_discretization_model_enabled
            dt = dt/(2.0*mesh.solid_1d_submesh.no_cells)
        end
    else
        dt = 0.48 * dt
    end

    if dt < 1.0E-6
        dt = 1.0E-6
        println("dt1: ", dts[1], "| dt2: ", dts[2], "| dt3: ", dts[3], "| dt4: "
                , dts[4])
        @warn "dt smaller dt_min, maybe something wrong here ..."

        if dt < 0
            error("negative timestep calculated, check bounding of pX functions")
        end
    end

    return dt
end

