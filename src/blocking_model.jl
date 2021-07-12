#=
Functions to estimate the blocking (factor) of scrap parts:

Here we assume a main line with length A_sum (which represents the cross 
sectional surface of the scrap shaft) and n_p lines with length A_p which 
represents the exposed surface of the individual scrap parts. 
So A_sum/(A_p * n_p) represents the occupied (area) fraction (occ_frac), 
determined by the radiation relevant area of earch scrap layer cell.

We distibute the n_p lines on the main line regulary and then randomly shift 
the offsets of the line in the range between half the gap size, between the 
original n_p lines. In this way we create a random distances, but avoid 
random clustering effects and very big gaps, which  propably will not occur 
in a real scrap bulk.

We then compare two randomly shiftet main lines, to each other. Under the 
assumption that they are infinitely close to each other we can then calculate
the overlap factor. From this we calculate a "seen amout" representing the amout
of surface, which contributes to radiation heat exchange between a 1d scrap cell 
i and its k neighboars (i-k:i-1,i+1:i+k).

As we see here due to symmetry between the lines (as each line has its own 
"movement space") the number of lines does not influence the results.

Some results are shown under ./ex/blocking_eval.jl
=#
################################################################################
################################################################################


function oc_blocking_over_layers(oc_fracs, n_layers::Int = 10)
    nps = [1] # only 1 needed - since symmetry ...
    mean_overlap_oc_frac = get_meanoverlap_over_oc_frac(nps, oc_fracs)
    layers = collect(0:1:n_layers)
    layers_between = layers .-1

    blocking_oc_over_layers = zeros(length(oc_fracs), length(layers))

    for (i,oc_frac) in enumerate(oc_fracs)
        for j = 1:length(layers)
            if j == 1
                blocking_oc_over_layers[i,j] = 1
            elseif j == 2
                blocking_oc_over_layers[i,j] = mean_overlap_oc_frac[i,end]
            else
                blocking_oc_over_layers[i,j] = (1-mean_overlap_oc_frac[i,end]).^layers_between[j] * mean_overlap_oc_frac[i,end]
            end
        end
    end

    return blocking_oc_over_layers
end


function get_meanoverlap_over_oc_frac(nps, oc_frac)
    A_sum = 1000.0     # fixed values
    ### investigate occupied_frac
    progress = Progress(length(nps)*length(oc_frac), 1) 

    update!(progress,0)
    mean_overlap_oc_frac = zeros(length(oc_frac),length(nps))

    for (j,n_p) in enumerate(nps)
        println("n_p = ", n_p)

        for i = 1:length(oc_frac)
            A_p_sum = A_sum * oc_frac[i]
            A_p = A_p_sum/n_p 
            mean_overlap_oc_frac[i,j] = calc_mean_overlap(A_sum, A_p, n_p, 20_000)
        end

    end

    return mean_overlap_oc_frac
end


function plot_overlap_example(A_sum, A_p, n_p)
    A_p_rel = A_p/A_sum

    p1 = get_rand_pos(1, A_p_rel, n_p)
    p2 = get_rand_pos(1, A_p_rel, n_p)

    lines_p1 = line_from_pos.(p1,A_p_rel)
    y_lines_p1 = y_line_from_pos(p1, 1.0)

    lines_p2 = line_from_pos.(p2,A_p_rel)
    y_lines_p2 = y_line_from_pos(p2, 2.0)

    lines_overlap, overlap = calc_overlap_plot(p1,p2,A_p_rel)
    y_lines_overlap = y_overlapline(lines_overlap, 1.5)

    lines_overlap = filter_missing(lines_overlap)
    y_lines_overlap = filter_missing(y_lines_overlap)

    x1_line_overlap, x2_line_overlap, y_line_overlap = 
            vertical_overlap_lines(lines_overlap, 1.0, 2.0)

    pl = plot(lines_p1, y_lines_p1, label ="", linecolor = :blue, linewidth = 3.0)
    plot!(pl, lines_p2, y_lines_p2,label ="", linecolor = :red, linewidth = 3.0)
    plot!(pl, lines_overlap, y_lines_overlap,label ="", linecolor = :green, 
            linewidth = 3.0)
    plot!(pl, x1_line_overlap, y_line_overlap,label ="", linecolor = :lightgray, 
            linestyle = :dash,linewidth = 1.0)
    plot!(pl, x2_line_overlap, y_line_overlap,label ="", linecolor = :lightgray, 
            linestyle = :dash,linewidth = 1.0)
    plot!(pl, xlim=(0,1), ylim = (0.9,2.1))
    plot!(pl, yticks=([1, 1.5, 2], [L"A_{rad,\, layer \, 1}", L"overlap", 
            L"A_{rad,\, layer \, 2}"]), bottom_margin = 5px, xlabel = "relative shaft cross section surface area")
    # display(pl)

    mean_mean_overlap = mean(overlap)

    println("Mean overlap: ", mean_mean_overlap)
    return pl
end

function filter_missing(lines::Vector{Vector{Union{Missing, Float64}}})

    lines_new = Vector{Vector{Float64}}(undef,0)

    for i = 1:length(lines)
        if !ismissing(lines[i][1]) && !ismissing(lines[i][2])
            push!(lines_new, [lines[i][1], lines[i][2]])
        end
    end

    return lines_new
end


function vertical_overlap_lines(lines_overlap, yp1, yp2)

    n = length(lines_overlap)
    x1_line_overlap = fill(Union{Float64,Missing}[0, 0],n)
    y_line_overlap = fill(Union{Float64,Missing}[0, 0],n)
    x2_line_overlap = fill(Union{Float64,Missing}[0, 0],n)

    for i = 1:length(lines_overlap)
        x1_line_overlap[i] = [lines_overlap[i][1], lines_overlap[i][1]]
        x2_line_overlap[i] = [lines_overlap[i][2], lines_overlap[i][2]]
        y_line_overlap[i] = [yp1, yp2]
    end

    return x1_line_overlap, x2_line_overlap, y_line_overlap
end

function line_from_pos(p_pos, A_p)
    return [p_pos - 0.5*A_p, p_pos + 0.5*A_p]
end

function y_line_from_pos(p_pos, yval)
    y_line_vec = ones(length(p_pos)) .* yval

    return y_val_vec.(y_line_vec)
end


function y_overlapline(overlap_lines, yval)
    n = length(overlap_lines)
    y_line_vec  = fill(Vector{Union{Missing,Float64}}(missing, 2),n)

    for i = 1:n
        if ismissing(overlap_lines[i][1])
            y_line_vec[i] = [missing, missing]
        else
            y_line_vec[i] = [yval, yval]
        end
    end

    return y_line_vec
end

function y_val_vec(yval)
    return [yval, yval]
end

function calc_overlap_plot(p1::Vector{T}, p2::Vector{T}, A_p::T) where T <: Real
    n = length(p1)
    lines_overlap = fill(Vector{Union{Missing,Float64}}(missing, 2),n)
    overlap = zeros(n)

    for i=1:n
        lines_overlap[i], overlap[i] = calc_overlap_plot(p1[i],p2[i],A_p)
    end

    return lines_overlap, overlap
end


function calc_overlap_plot(p1::Real, p2::Real, A_p::Real)
    overlap_tol = 1E-8

    plow, overlap = range_p_overlap(range_p(p1,A_p), range_p(p2,A_p), A_p)
   
    if !ismissing(plow)
        overlap_line = [plow, plow+overlap]
    else
        overlap_line = [missing, missing]
    end

    overlap = overlap / A_p  # in percent

    return overlap_line, overlap
end

function calc_mean_overlap(A_sum::Float64, A_p::Float64, n_p::Int64, n_runs = 20)
    mean_mean_overlap = zeros(n_runs)

    @inbounds for i = 1:n_runs
        p1 = get_rand_pos(A_sum, A_p, n_p)
        p2 = get_rand_pos(A_sum, A_p, n_p)
        overlap = calc_overlap(p1,p2,n_p,A_p)
        mean_mean_overlap[i] = mean(overlap)
    end

    mean_overlap = mean(mean_mean_overlap)
    return mean_overlap
end


function calc_overlap(p1, p2, n_p, A_p)
    overlap_tol = 1E-8
    overlap = zeros(Float64, n_p)

    @inbounds for i = 1:n_p
        local plow
        plow, overlap[i] = range_p_overlap(range_p(p1[i],A_p), range_p(p2[i],A_p), A_p)
        # in percent
        overlap[i] = overlap[i] / A_p
    end

    if ( (sum(overlap .< -overlap_tol) >= 1) ||  (sum(overlap .> 1+overlap_tol) >= 1) )
        println( overlap[(overlap .< 0) .| (overlap .> 1)] )
        error("Error in overlap calulation")
    end

    return overlap
end

function range_p(p_pos, A_p)
    return (p_pos - 0.5*A_p, p_pos + 0.5*A_p)
end

function range_p_overlap(p1_range, p2_range, A_p)
    if p1_range[1] > p2_range[1]
        # p1_min > p2_min
        if p1_range[1] > p2_range[2]
            # p1_min > p2_max
            return missing, 0
        end

        return p1_range[1], (A_p - (p1_range[1] - p2_range[1]))

    elseif p1_range[1] < p2_range[1]
        # p2_min > p1_min
        if p1_range[2] < p2_range[1]
            # p1_min > p2_max
            return missing, 0
        end

        return p2_range[1],(A_p - (p2_range[1] - p1_range[1]))
    else
        # equal
        return p1_range[1], A_p
    end
end

function get_rand_pos(A_sum, A_p, n_p)
    off_min = 0.5 * A_p
    off_gap =  (A_sum-A_p*n_p)/n_p
    off_regular = off_gap+2*off_min
    p_pos = collect(off_min+off_gap/2:off_regular:A_sum-off_gap/2)
    dist = zeros(Float64, n_p)

    @inbounds for i=1:n_p
        p_pos[i] = p_pos[i] + rand(-off_gap/2:(off_gap/50):off_gap/2)
        
        if i>1 
            dist[i] = p_pos[i] - p_pos[i-1] 
        end
    end

    error_check = falses(n_p)
    error_check = dist .< off_min-1E-16

    if sum(error_check[2:end]) != 0
        error("Something went wrong generating random pos")
    end

    return p_pos
end

