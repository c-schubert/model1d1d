module model1d1d
using Base: AbstractFloat
using Dates
using QuadGK
using LinearAlgebra
using ProgressMeter
using Statistics
using LaTeXStrings
# suppressses the GKS QtTerm Windows from appreaing when using plots for savefig
ENV["GKSwstype"] = "nul" 
using Plots
using Plots.PlotMeasures
gr()
using Printf
using Colors
using StaticArrays

const Tmax = 2000.0
const Twarn = 1900.0
const Trange = LinRange(0.0,Tmax+1, 10_000)  
const sigma = 5.670374419E-8
const dt_max = 0.1

# use for logging inside module files
# initialize using new_log to custom logfile location (name)
m1d1d_logfile = "m1d1dlogfile.log"

include("print_log.jl")
include("front_radiation_types.jl")
include("matprop_types.jl")
include("matprop_functions.jl")
include("heat_transfer_types.jl")
include("heat_transfer_functions.jl")
include("model1d1d_types.jl")
include("mesh_types.jl")
include("varcalc.jl")
include("faster_vec_ops.jl")
# solver function using self written explicit euler methods
# may be instable (should be used for testing only)
# include("solver_1d1d_ee.jl")
# Main solver function using DifferentialEquations.jl
include("solver_1d1d_diffeq.jl")
include("blocking_model.jl")
include("post_process.jl")



export print_log, new_log, test_global_logfilename
export solve_1d1d_diffeq
export SimSettings

end