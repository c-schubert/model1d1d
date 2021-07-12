"""
Functions to return MaterialPropertyOfX for X appropriate to type 
MaterialPropertyOfX
"""

"""
interp1_bounded:
bounded linear interpolation of values xarr, yarr for value v
"""
function interp1_bounded(xarr::Vector{<:Real},yarr::Vector{<:Real} , v::Real)
    #xarr musst be sorted low to hight
    tol = 1E-8
    res = 0.0

    if v <= xarr[1]
        res =  yarr[1]
    elseif v >= xarr[end]
        res =  yarr[end]
    else   
        xlow = findfirst(f -> f .> (v+tol), xarr)

        res = ((yarr[xlow] - yarr[xlow-1])/(xarr[xlow] - xarr[xlow-1]) 
            * (v-xarr[xlow-1]) + yarr[xlow-1])
    end

    return res
end

"""
interp1_extra:
linear extrapolated linear interpolation of values xarr, yarr for value v
"""
function interp1_extra(xarr::Vector{<:Real},yarr::Vector{<:Real} , v::Real)
    #xarr musst be sorted low to hight
    tol = 1E-8
    
    if v <= xarr[1]
        xlow = 2
        res = ((yarr[xlow] - yarr[xlow-1])/(xarr[xlow] - xarr[xlow-1]) 
        * (v-xarr[xlow-1]) + yarr[xlow-1])

    elseif v >= xarr[end]
        xlow = length(xarr)
        res = ((yarr[xlow] - yarr[xlow-1])/(xarr[xlow] - xarr[xlow-1]) 
        * (v-xarr[xlow-1]) + yarr[xlow-1])
    else
        xlow = findfirst(f -> f .> (v+tol), xarr)

        res = ((yarr[xlow] - yarr[xlow-1])/(xarr[xlow] - xarr[xlow-1]) 
            * (v-xarr[xlow-1]) + yarr[xlow-1])
    end

    return res
end

"""
In the follwowing the implementations for the different MaterialPropertyOfX are given:
"""
# import Base

Base.iterate(x::MaterialPropertyOfX) = (x, nothing)
Base.iterate(x::MaterialPropertyOfX, ::Any) = nothing
Base.isempty(x::MaterialPropertyOfX) = false
Base.length(x::MaterialPropertyOfX) = 1


"""
Sepcific vector versions of functions allow for little optimization for some types
TODO: Add Variantes 
"""
function pX(p::T1, X::T2)::T2 where T2 <: Vector{<:Real} where T1 <: MaterialPropertyOfX
    return pX.(p, X)
end

function pX_bounded(p::T1, X::T2)::T2 where T2 <: Vector{<:Real} where T1 <: MaterialPropertyOfX
    return pX.(p, X)
end

function pX_bounded(p::T1, X::T2)::T2 where T2 <: Vector{<:Real} where T1 <: BoundedMaterialPropertyOfX  
    return  pX_bounded.(p, X)
end

function pX(p::T1, X::T2)::T2 where T2 <: Vector{<:Real} where T1 <: ConstantMaterialPropertyOfX
    return fill( pX(p, X[1]), length(X))
end

function pX_bounded(p::T1, X::T2)::T2 where T2 <: Vector{<:Real} where T1 <: ConstantMaterialPropertyOfX
    return pX(p, X)
end


"""
MaterialPropertyOfX types which have no bounds property should be able 
to bound anyhow
"""
function pX_bounded(p::T1, X::T2)::T2 where T2 <: Real where T1 <:MaterialPropertyOfX
    p_val = pX(p, X)
    
    return p_val
end

function pX_bounded(p::T1, X::T2)::T2 where T2 <: Real where T1 <:BoundedMaterialPropertyOfX
    p_val = pX(p,X)

    if !isnothing(p.bounds)
        if p_val < p.bounds[1]
            p_val = p.bounds[1]
        end

        if p_val > p.bounds[2]
            p_val = p.bounds[2]
        end
    end

    return p_val
end

function pX(p::T1, X::T2)::T2 where T2 <: Real where T1 <:MaterialPropertyOfX
    error("This function should not be called ...")
    return X
end

function pX(p::ConstantMaterialPropertyOfX{T}, X::T)::T where T <: Real
    return p.val
end

function pX(p::TabulatedMaterialPropertyOfX, X::T)::T where T <: Real
    index = floor(Int64, ((X - p.range.start) + p.rangestep + 1E-16) / p.rangestep)

    return @boundscheck p.val_arr[index]
end

function pX(p::LinearInterpMaterialPropertyOfX{T}, X::T)::T where T <: Real
    if p.bounded
        return interp1_bounded.(p.x,p.val, X)
    else
        return interp1_extra.(p.x,p.val, X)
    end
end

function pX(p::LinearMaterialPropertyOfX{T}, X::T)::T where T <: Real
    return p.m * X + p.y1
end

function pX(p::PolynomialMaterialPropertyOfX{T1, T2}, X::T1)::T1 where {T1 <: Real, T2 <: Integer}
    p_val = 0.0
    for i = 1:1:p.n
        p_val += p.a[i] * X^(i-1)
    end
    return p_val
end

function pX(p::ReciprocalMaterialPropertyOfX{T}, X::T)::T where T <: Real
    return p.m * 1.0 /X + p.y1
end

function pX(p::ExponentialMaterialPropertyOfX{T}, X::T)::T where T <: Real
    return p.a * exp(p.b * x) + p.y1
end

function pX(p::LogarithmicMaterialPropertyOfX{T}, X::T)::T where T <: Real
    return p.a * ln.(X) + p.y1
end

"""
num_integrate_pX:
in general use numerical integration
(uses QuadGK.jl)
"""
function num_integrate_pX(p::MaterialPropertyOfX, 
        X_from::T, X_to::T)::T where T <: Real 

    H, err = quadgk.( x -> pX(p, x), X_from, X_to)

    if err > 1e-3
        error("Error integration of x is high")
    end
    return H
end

"""
integrate_pX(p::PolynomialMaterialPropertyOfX, X_from::Real, 
        X_to::Real):

Analytical integration of polynomial material properties
"""
function integrate_pX(p::PolynomialMaterialPropertyOfX{T, <:Integer}, X_from::T, 
        X_to::T)::T where T <: Real 
        
    local  p_val_from::T 
    local p_val_to::T
    p_val_from = 0.0
    p_val_to = 0.0

    for i = 1:1:p.n
        p_val_from += 1.0/(i) * p.a[i] * X_from^(i)
        p_val_to += 1.0/(i) * p.a[i] * X_to^(i)
    end

    return (p_val_to - p_val_from)
end


"""
integrate_pX(p::ConstantMaterialPropertyOfX, X_from::Real, 
        X_to::Real):

Analytical integration of constant material properties
"""
function integrate_pX(p::ConstantMaterialPropertyOfX{T}, X_from::T, 
        X_to::T)::T where T <: Real 
        
    return (p.val * (X_to - X_from))
end


"""
integrate_pX(p::TabulatedMaterialPropertyOfX, X_from::Real, 
        X_to::Real):

Analytical integration of Tabulated material properties
"""
function integrate_pX(p::TabulatedMaterialPropertyOfX, X_from::T, 
        X_to::T)::T where T <: Real 

    ifrom = floor(Int64, (X_from - (p.range.start - p.rangestep)) / p.rangestep)
    ito = floor(Int64, (X_to - (p.range.start - p.rangestep)) / p.rangestep)

    int = 0.0
    if ito > ifrom
        @boundscheck for i = ifrom:1:ito
            int += p.val_arr[i] * p.rangestep
        end
    else
        @boundscheck for i = ito:-1:ifrom
            int -= p.val_arr[i] * p.rangestep
        end
    end

    return int
end


"""
integrate_pX(p::MaterialPropertyOfX, X_from::Real, 
        X_to::Real):

Fall back to numerical integration until specific version is available!

TODO: Impelement integration functions for other MaterialPropertyOfX
"""
function integrate_pX(p::MaterialPropertyOfX, X_from::T, 
    X_to::T)::T where T <: Real 

    return num_integrate_pX(p, X_from, X_to)
end

