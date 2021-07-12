"""
Various functions for different calculations
"""

"""
surface_quader:
Calculates the surface of a quader.
"""
function surface_quader(W::Real, L::Real, T::Real)

    return (2*L*W+2*T*L+2*T*W)
end


"""
volume_quader:
Calculates the surface of a quader.
"""
function volume_quader(W::Real, L::Real, T::Real)

    return (W*L*T)
end


"""
ct_and_cA_quader:

Returns the characteristic thickness of a quader.
The characteristic thickness is the greates dimension normal to the greates 
surface plane (characteristic cross sectional area see bellow) of the quader
"""
function char_thickness_quader(dims::Vector{T}) where T <: Real

    if length(dims) != 3
        error("char_thickness_quader expects an array length 3")
    end

    A = [dims[1]*dims[3], dims[1]*dims[2], dims[2]*dims[3]]
    n = [dims[2], dims[3], dims[1]]

    return maximum(n[A .== maximum(A)])
end


"""
char_cross_section_area_quader:

Returns the characteristic cross sectional area of a quader.
"""
function char_cross_section_area_quader(dims::Vector{T}) where T <: Real

    if length(dims) != 3
        error("char_cross_section_area_quader expects an array length 3")
    end

    A = [dims[1]*dims[3], dims[1]*dims[2], dims[2]*dims[3]]
    n = [dims[2], dims[3], dims[1]]

    return maximum(A)
end


"""
 surface_area_without_contact():

Calculates the maximal free surface area of pieces with surface 
"surface_of_piece" for a given bulk density "bulk_density"
"""
function surface_area_without_contact(bulk_density::Real, weight_of_piece::Real, surface_of_piece::Real )

    no_pieces = (bulk_density/scrap_weight_of_piece)
    return (no_pieces * surface_of_piece)
end


"""
hf_radiation_inf(eps::Real, A::Real,T1::Real, T2::Real)

returns the heat flux from T1 to T2 and this formula is valid if A2 is much 
bigger than A1 if epsilon2 >> 0.
"""
function hf_radiation_inf(eps::Real, A1::Real, T1::Real, T2::Real)

    return (eps * A1 * sigma * (T1^4 - T2^4))
end

"""
hf_radiation_par_surf(eps1::Real, eps2::Real, T1::Real, T2::Real)

returns the heat flux (W/m² K) between to infinite parallel surfaces from 1 to 2
"""
function hf_radiation_par_surf(eps1::Real, eps2::Real, T1::Real, T2::Real)

    return (sigma * (T1^4 - T2^4) / (1/eps1 + 1/eps2 - 1))
end


"""
hf_radiation_par_surf_b_weighted(eps::Real, b::Real, T1::Real, T2::Real)

returns the heat flux (W/m² K) between to equally sized parallel surfaces from 1 
to 2 which is weighet by a factor of b
"""
function hf_radiation_par_surf_b_weighted(eps1::Real, eps2::Real, b::Real, T1::Real, T2::Real)

    return ( b * hf_radiation_par_surf(eps1, eps2, T1, T2) )
end


function H_of_T(cp_true::MaterialPropertyOfX, T::Real)::Real

    if T > Tmax - eps(typeof(T)) || T < 0 - eps(typeof(T))   
        error("value " * string(T), " out of bounds for H_of_T!")
    end

    return integrate_pX(cp_true, 0.0, T)
end

"""
T_of_H:
Caluculates the temperature T out of the enthalpy curve H(T), using bin search 
so thats it should be relativly fast ...

If we are using T(H) curve to estimate T is it so that H has to have no 
minima or maxima? - it should, as cp is deterimed by heating curves ...
"""
function T_of_H(H_T::Vector{<:Real}, T_H::Vector{<:Real}, H::Real)

    i_low = 1
    i_high =  length(H_T)
    while true
        i = ceil(Int64, (i_low + i_high - 1)/2)

        if H < H_T[i]
            i_high = i
        else
            i_low = i
        end

        if (i_high - i_low) == 1
            break;
        end
    end

    return  ((T_H[i_high] - T_H[i_low])/(H_T[i_high] - H_T[i_low]) 
                    *(H - H_T[i_low]) + T_H[i_low])
end

function T_of_H(H_T::Vector{<:Real}, T_H::Vector{<:Real}, H::Vector{<:Real})::Vector{<:Real}

    T = zeros(length(H))

    for i = 1:length(H)
        T[i] =  T_of_H(H_T, T_H, H[i])
    end

    return T
end


"""
check_monolithic_growing(p::Vector{<:Real})::Nothing
Checks if values in Vector p are monolithicly growing, used for checking 
enthalpy curves, see above.
"""
function check_monolithic_growing(p::Vector{<:Real})::Nothing

    for i = 2:1:length(p) 
        ((p[i] < p[i-1]) ? warning("Non growing values detected for H") : 
                continue)
    end
end


"""
harm_mean(x::Real...)::Real
Simple (naive) implementation of harmonic mean average
"""
function harm_mean(x::Real...)::Real
    n = length(x)
    denom = 0
    for xi in x
        denom += 1/xi
    end

    return (n/denom)
end


"""
harm_mean(x1::T,x2::T)
optimization if harm mean only gets to values - does not need any 
memory allocation ...
"""
function harm_mean(x1::T,x2::T)::T where T <: Real
    denom = 1/x1 + 1/x2

    return (2/denom)
end


"""
print_struct(T::Type, strct)
Function to print a struct and its properties
"""
function print_struct(T::Type, strct)
    for fn in fieldnames(T)
        fval = getfield(strct, fn)

        println("VAR:\n", string(fn), "\nTYPE:\n{", typeof(fval),"}\nVAL:\n", fval,"\n\n")
    end
end
