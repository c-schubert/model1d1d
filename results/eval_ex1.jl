ENV["GKSwstype"] = "nul"
using model1d1d
# using Plots
# using Plots.PlotMeasures
using JLD2

function eval_all_times(;folder::AbstractString=pwd(), img_folder::AbstractString=joinpath(pwd(), "img"))

    folder = joinpath(pwd(),"results")
    casenames = String[]

    files = readdir(folder)
    jld2_files = files[endswith.(files, ".jld2")]

    @show jld2_files

    for jld2_file in jld2_files

        casename = jld2_file[1:end-5]

        if !isdir(joinpath(img_folder, casename))
            println("Generating images for ", casename)

            @load joinpath(folder, jld2_file) model desc results 

            println("Casename: ", casename, "\n")

            model1d1d.print_case_info(model)

            model1d1d.post_process_results(@__DIR__, casename, model, results, 
                    ylim_cust = (0.0,1400.0))
        else
            println("Result image folder allready exsits, skipping!")
        end
    end
end

eval_all_times()

