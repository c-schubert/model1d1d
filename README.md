# 1d1d Scrap Preheating Model

This model was developed to modell the preheating of scrap in an semi-continuous 
EAF scrap charging shaft.


## Setup and Usage

> **Only tested with julia 1.5.3 and 1.6.1**

> _All files should be run from within the project folder!_

> `./src` contains the essential code of the model

### 0. Download and install julia
  1. Download from: www.julialang.org
  2. Add the folder containing the julia.exe to system variable: `Path` 

### 1. Setup Environment & Run Scripts
  1. Open a Terminal 
  2. Goto the folder of this project containing the `Project.toml` file.
  3. Open Julia
  4. Switch to package mode with  `]`
  5. In pkg mode activate current environment with `activate .`
  6. Use `instantiate` to install packages defined in the `Project.toml` file (first time only)
  7. Switch back to julia REPL mode (backspace)
  8. run desired script from within julia for example with `include("./ex/eval_ex.1")`


# Future Development Ideas (TODO List ...)
This is a pretty long list btw ...

> **Discalimer:** This project is currently not in active development, as there 
> is no funding to work on this topic at the moment. So you can read the list 
> more as a what could be done, instead of what will be done. If you are interessted
> if future developments of this model feel free to contact us over at www.iob.rwth-aachen.de

## Performance Optimization 

#### Spliting 1d-1d Modelling
We could split the time stepping between the heat conduction problem in the scrap pieces, 
and the heat transport in the overall shaft.

Therefore we could use dt_shaft and dt_scrap (probably dt_shaft will be much higher) 
and only adapt boundaries for heat conduction in the scrap after each dt_shaft. 

Futhermore this could have the benefits that we could decouple the methods used to 
solve the problems. 

This may lead to some inaccuracies if dt_scrap grows to large, so we have to be 
careful here.

## Usability
Extend Login to make cases more reproducible - maybe input file format?

### Settings/Property Files
Offload the settings which are right not done in the `prepare_data.jl` to i.e. 
specific json files.

## Features
Possible features that could be implemented.

### More specific front side radiation modell
In a previous version of this model (written in Matlab) we had an surface2-surface 
(s2s) model, but it was very specific for a certain problem. Generalisation of 
such a model is hard - maybe we could interface to FARADS3D at some point ...

### Scrap Piece Conduction Modell
Right now we ignoring the heat conduction between the individual scrap pieces. 
Maybe we can come up with a simple model based on contact areas and contact heat 
resistances. But this should not be to relevant in the most use cases ...

### Convection Models (HTC)
Incorporate some models to relate changes of the off-gas flow rate to changes of 
the heat transfer coefficients (htc) as it depends on the Reynolds-Number (Re). 
Therefore we should considering evaluating changes of the off-gas properties 
(Pr, lam, nu) to the htc, as they influence Re and Nu (Nusselt-Number).

### Enhance Scrap Piece Radiation Modelling
The modeling right now is very simple and uses many assumptions, can we use other 
models like FARADS2D to enhance the radiation modelling a bit without getting 
to complex?

### More Integration Functions
Add more specific integration function for all `MaterialPropertyOfX`s in 
`matprop_functions.jl`

### Scrap Carbon Modell 
Implement some type of carbon exchange model for bulk and off-gas to also use 
the model for post combustion optimization calculations.

### Auto Optimization
Allowing data based auto optimization (maybe NN or simple optimization methods) 
of modelling coefficients.

### Make Model 2-3 dimensional
Extend 1d1d model to 2d1d or 3d1d model, so that spacially different scrap 
compositions and properties can be used.

## General Cleanup
  - make solve1d1d(m::M1d1d) more modular
  - can we get rid of manual Float64 type annotations in solver1d1d?
  - enhance precompile time

## Documentation
More and more better!

