
function unsafe_update_v1_to_v2!(v1::Vector{T}, v2::Vector{T}, n::Int) where T<:Real
    @inbounds @simd for i = 1:n
        v1[i] = v2[i]
    end
end

function unsafe_update_vv1_to_vv2!(vv1::Vector{Vector{T}}, vv2::Vector{Vector{T}}, n1::Int , n2::Int) where T<:Real

    @inbounds for i = 1:n1
        @simd for j = 1:n2
            vv1[i][j] = vv2[i][j]
        end
    end
end

function my_copy(vec::Vector{T})::Vector{T} where T<:Real
    n = length(vec)
    new_vec = Vector{T}(undef, n)

    @inbounds @simd for i = 1:n
        new_vec[i] = vec[i]
    end

    return new_vec
end

# test
# a = [[rand(1:5) for j in 1:5] for i in 1:200]
function my_deep_copy_unsafe(vec_vec::Vector{Vector{T}}, n1::Int , n2::Int)::Vector{Vector{T}} where T<:Real
    new_vec = Vector{Vector{T}}(undef, n1)

    @inbounds for i = 1:n1
        new_vec[i] = Vector{T}(undef, n2)

        @simd for j = 1:n2
            new_vec[i][j] = vec_vec[i][j]
        end
    end

    return new_vec
end
