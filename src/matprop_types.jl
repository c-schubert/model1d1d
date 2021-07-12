"""
Implementation of types to describe the behavior of a material from a dependent 
variable.
"""

#= #############################################################################
Material Properties "Law" Types
############################################################################# =#

"""
MaterialPropertyOfX:
describes a whatever (on x) depended material property (p)

We use several subtypes to allow for more accurate material property fitting
"""
abstract type MaterialPropertyOfX end;
abstract type BoundedMaterialPropertyOfX <: MaterialPropertyOfX end;


"""
ConstantMaterialPropertyOfX:
property that has not dependency on X
"""
struct ConstantMaterialPropertyOfX{T <: AbstractFloat} <: MaterialPropertyOfX
    name::String
    unit::String
    val::T

    ConstantMaterialPropertyOfX(val::T;
                                name::String = "empty",
                                unit::String = "empty"
                                ) where T <: AbstractFloat = 
                                new{T}(name,unit,val)
end


"""
TabulatedMaterialPropertyOfX:
tabulated material property

Using tabulated material properties does not make many sense right now as the 
computation of material properties itself is very cheap (even polymial)
"""
struct TabulatedMaterialPropertyOfX{T1 <: LinRange{<:AbstractFloat}, 
    T2 <: AbstractFloat, T3 <: Vector{<:AbstractFloat}} <: MaterialPropertyOfX
    name::String
    unit::String
    range::T1 # .len .stop .start
    rangestep::T2
    val_arr::T3

    TabulatedMaterialPropertyOfX(range::T1, 
                                val_arr::T2;
                                name::String = "empty",
                                unit::String = "empty"
                                ) where {T1 <: LinRange{<:AbstractFloat},
                                T2 <: Vector{<:AbstractFloat}} = 
                                new{T1, typeof(((range.stop-range.start)/range.len)),T2}( 
                                name, unit, range, 
                                ((range.stop - range.start) / range.len), 
                                val_arr)
end


function TabulatedMaterialPropertyOfX(MP::MaterialPropertyOfX, 
        range::LinRange{<:AbstractFloat}, bounded::Bool = true)
    # we get T of MaterialPropertyOfX
    # MPtype = typeof(MP)
    # MPdatafieldtype = MPtype.parameters[1]
    # val_arr = zeros(MPdatafieldtype, range.len)

    if bounded
        val_arr= pX_bounded.(MP, range)
    else
        val_arr = pX.(MP, range)
    end

    return TabulatedMaterialPropertyOfX(range, val_arr, name = MP.name, unit = MP.unit)
end

"""
LinearInterpMaterialPropertyOfX:
Implements type for linear interpolated materialproperties. Values are lineary 
interpolated between die nearest values. If bounded is true no linear 
extrapolation will be done and the last value will be used as constant 
extrapolation.
"""
struct LinearInterpMaterialPropertyOfX{T <: AbstractFloat} <: MaterialPropertyOfX
    name::String
    unit::String
    bounded::Bool
    x::Vector{T}
    val::Vector{T}

    LinearInterpMaterialPropertyOfX(x::Vector{T}, val::Vector{T};
                                    bounded::Bool = false,
                                    name::String = "empty",
                                    unit::String = "empty"
                                    ) where T <: AbstractFloat = 
                                    new{T}(name,unit,bounded,x,val)
end

"""
LinearMaterialPropertyOfX:
Implements material property based on linear law

p(x) = m*x + y1
"""
struct LinearMaterialPropertyOfX{T <: AbstractFloat} <: BoundedMaterialPropertyOfX
    name::String
    unit::String
    m::T
    y1::T
    bounds::Union{Nothing, Tuple{T,T}}

    LinearMaterialPropertyOfX(m::T, y1::T;
                            name::String = "empty",
                            unit::String = "empty", 
                            bounds::Union{Nothing, Tuple{T,T}} = nothing
                            ) where T <: AbstractFloat = 
                            new{T}(name,unit,m,y1, bounds)
end

"""
PolynomialMaterialPropertyOfX:
Implements material property based on polynomial law from till degree n

p(x) = a[1] + a[2]*x + a[3] * x^2 + ... + a[n] * x^n
"""
struct PolynomialMaterialPropertyOfX{T1 <: AbstractFloat, T2 <: Integer}  <: BoundedMaterialPropertyOfX
    name::String
    unit::String
    n::T2
    a::Vector{T1}
    bounds::Union{Nothing, Tuple{T1,T1}}

    PolynomialMaterialPropertyOfX(a::Vector{T};
                                  name::String = "empty",
                                  unit::String = "empty",
                                  bounds::Union{Nothing, Tuple{T,T}} = nothing
                                ) where T <: AbstractFloat = 
                                new{T, typeof(length(a))}(name,unit,length(a),a,bounds)

end

"""
ReciprocalMaterialPropertyOfX:
Reciprocal material property

p(x) = m * 1/x + y1
"""
struct ReciprocalMaterialPropertyOfX{T <: AbstractFloat} <: BoundedMaterialPropertyOfX
    name::String
    unit::String
    m::T
    y1::T
    bounds::Union{Nothing, Tuple{T,T}}

    ReciprocalMaterialPropertyOfX(m::T, y1::T;
                                    name::String = "empty",
                                    unit::String = "empty",
                                    bounds::Union{Nothing, Tuple{T,T}} = nothing
                                    ) where T <: AbstractFloat = 
                                    new{T}(name,unit,m,y1,bounds)
end

"""
ExponentialMaterialPropertyOfX:
Exponential material property

p(x) = a * exp(b*x) + y1
"""
struct ExponentialMaterialPropertyOfX{T <: AbstractFloat} <: BoundedMaterialPropertyOfX
    name::String
    unit::String
    a::T
    b::T
    y1::T
    bounds::Union{Nothing, Tuple{T,T}}

    ExponentialMaterialPropertyOfX(a::T,b::T, y1::T;
                                    name::String = "empty",
                                    unit::String = "empty",
                                    bounds::Union{Nothing, Tuple{T,T}} = nothing
                                    ) where T <: AbstractFloat = 
                                    new{T}(name,unit,a,b,y1,bounds)
end

"""
LogarithmicMaterialPropertyOfX:
Logarithmic material property

p(x) = a * ln(x) + y1
"""
struct LogarithmicMaterialPropertyOfX{T <: AbstractFloat} <: BoundedMaterialPropertyOfX
    name::String
    unit::String
    a::T
    y1::T
    bounds::Union{Nothing, Tuple{T,T}}
    
    LogarithmicMaterialPropertyOfX(a::T, y1::T;
                                    name::String = "empty",
                                    unit::String = "empty",
                                    bounds::Union{Nothing, Tuple{T,T}} = nothing
                                    ) where T <: AbstractFloat = 
                                    new{T}(name,unit,a,y1,bounds)
end

#= #############################################################################
Final Material Types
############################################################################# =#

"""
IncompressibleFluid:
MaterialPropertyOfX here X should represent the temperature

"""
struct IncompressibleFluidThermalMaterial{T1 <: MaterialPropertyOfX,
        T2 <: MaterialPropertyOfX,
        T3 <: MaterialPropertyOfX}

    name::String
    rho::T1
    cp_true::T2
    lambda::T3
end

"""
SolidThermalMaterial:
MaterialPropertyOfX here X should represent the temperature

"""
struct SolidThermalMaterial{T1 <: MaterialPropertyOfX,
        T2 <: MaterialPropertyOfX, T3 <: MaterialPropertyOfX, 
        T4 <: MaterialPropertyOfX}

    name::String
    rho::T1
    cp_true::T2
    lambda::T3
    epsilon::T4
end
