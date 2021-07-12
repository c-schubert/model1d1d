function plot_res_at_ind(model::M1d1d, res::Vector{Results1d1d},
    t_ind::Integer, degC::Bool = false)

    ylabeltxt = "Temperature in K"
    TinX = T -> T
    if degC
        TinX = T -> (T .- 273.15)
        ylabeltxt = "Temperature in °C"  
    end

    t = res[t_ind].t

    x = pp_get_bulk_mesh_coords(res[t_ind], model)
    # println("Plotting results for time ", seconds_to_min_seconds(t), " seconds")

    scrap_conduction = false
    i = 1
    t_cum_phase_dur = 0.0
    while t_cum_phase_dur <= t && i <= length(model.sim_phase)
        t_cum_phase_dur += model.sim_phase[i].duration
        scrap_conduction = model.sim_phase[i].submodels.scrap_discretization_model_enabled
        i += 1
    end

    if scrap_conduction          
        pl = plot(x, TinX(res[t_ind].T.s), label = "T solid (mean)");
        T_sub_max = [res[t_ind].T.s_sub[i][1] for i in 1:length(res[end].T.s_sub)]
        T_sub_min = [res[t_ind].T.s_sub[i][end] for i in 1:length(res[end].T.s_sub)]
        plot!(pl, x, TinX(T_sub_max), label = "T solid (max)");
        plot!(pl, x, TinX(T_sub_min), label = "T solid (min)");
    else
        pl = plot(x, TinX(res[t_ind].T.s), label = "T solid");
    end

    plot!(pl, title = seconds_to_min_seconds(t),
            font = font(12,"Helvetica", :black), xlabel = "x in m", 
            ylabel = ylabeltxt)

    plot!(pl, x, TinX(res[t_ind].T.f), label = "T fluid");

    return pl;
end


function pp_get_bulk_mesh_coords(res::Results1d1d, model::M1d1d)
    t = res.t
    mesh = nothing
    dur = 0.0

    for pid = 1:length(model.sim_phase)
        dur += model.sim_phase[pid].duration

        if t <= dur
            mesh = M1d1dMesh(model, pid)
            break
        end
    end

    if isnothing(mesh)
        error("No mesh!")
    end

    return [i*mesh.solid_cells.len for i=0:1:mesh.no_cells-1]
end


function seconds_to_min_seconds(seconds::Real)
    min = floor(Int64, seconds/60)
    sec = round(Int64,mod(seconds,60))

    return string(min) * " min " * string(sec) * " sec"
end


function post_process_results(pp_path::String, casename::String, model::M1d1d,
        res::Vector{Results1d1d};  in_deg_C = true, ylim_cust = nothing, dpi_pic=100)

    pp_img_path = joinpath(pp_path,"img",casename)
    clear_check_folder(pp_img_path)
    println("Generating images for all saved timesteps, saving to folder ", pp_img_path)
    progress = Progress(length(res), 1) 
    for i =1:length(res)
        pl = model1d1d.plot_res_at_ind(model, res, i, in_deg_C);
        plot!(pl, dpi=dpi_pic, fmt=:png);
        isnothing(ylim_cust) ? nothing : plot!(pl, ylim = ylim_cust)
        png(pl, joinpath(pp_img_path, 
                casename * "-" *replace(seconds_to_min_seconds(res[i].t)," " => "_") * ".png") );
        next!(progress)
    end

    return
end


function print_case_info(m::M1d1d)
    println("-------------------------------------------")
    for i=1:1:length(m.bulk)
    
        println("Bulk number: ", i)
        println("Bulk mass flow: ", m.bulk[i].mass_flow)
        println("Bulk density: ", m.bulk[i].density)
        println("Bulk porosity: ", m.bulk[i].porosity)
        println("Bulk cross section area: ", m.bulk[i].cross_sect_surf_area)

        println("Scrap part surface area: ", m.bulk[i].scrap_props.surf_area)
        println("Scrap part thickness: ", m.bulk[i].scrap_props.char_thickness)
        println("Scrap part volume: ", m.bulk[i].scrap_props.vol)

        mesh = M1d1dMesh(m,i)
        printInfo(mesh)
        println("\n\n")
    end
    
end

    
function clear_check_folder(path::String)
    # check if folder is empty / exists
    if !ispath(path) && !isfile(path)
            mkpath(path)
    elseif isfile(path)
            error("cannot make path, because path is allready a file location")
    end
end