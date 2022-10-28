

mutable struct SafeInterpolator
    itp::Union{Interpolations.ScaledInterpolation}
    itp_safety::Union{Interpolations.ScaledInterpolation}
end

function interpolate(interpolator::SafeInterpolator, index, tol::Float64=0.5)
    default_interpolation = exp(Gr, origin, interpolator.itp(index))
    # find the interpolation points
    range = interpolator.itp.ranges[1]
    right = 0
    left = 0
    for i in range
        if i == index
            return default_interpolation # this is a sample point
        end
        if i > index 
            right = i
            break
        end
        left = i
    end
    safety_interpolation = exp(Gr_safety, origin_safety, interpolator.itp_safety(index))

    sim_right = abs(dot(default_interpolation[:,n_eigenfunctions], exp(Gr_safety, origin_safety, interpolator.itp_safety(right))[:,n_eigenfunctions+1]))
    sim_left = abs(dot(default_interpolation[:,n_eigenfunctions], exp(Gr_safety, origin_safety, interpolator.itp_safety(left))[:,n_eigenfunctions]))

    if sim_left > sim_right
        optimal_interpolation = hcat( default_interpolation[:,1:end-1], safety_interpolation[:,n_eigenfunctions])
    else
        optimal_interpolation = hcat( default_interpolation[:,1:end-1], safety_interpolation[:,n_eigenfunctions + 1])
    end
    return optimal_interpolation
end

function interpolate(interpolator::Union{Interpolations.ScaledInterpolation, Interpolations.Extrapolation, Interpolations.BSpline}, index)
    return interpolator(index)
end


function interpolate_complete( itp_sample_points, itp_sample_points_parameter, itp_parameter_points, origin, n_dim)
    Gr =  Grassmann(size(itp_sample_points)[1], n_dim)
    # 1. project manifold points onto the tangentspace
    V = Matrix{Float64}[]  # eigenvectors transformed to the tangentspace of the origin
    map(i -> push!(V, log(Gr, origin, i)), itp_sample_points)

    # 3. interpolate on the tangentspace
    itp = interpolate(V, BSpline(Linear()))
    sitp = scale(itp, itp_sample_points_parameter)
    interpolated_tangentspace = [sitp(parameter) for parameter in itp_parameter_points]

    # 3. retract from the tangentspace back to the manifold
    interpolated_solution = map(i -> exp(Gr, origin, i), interpolated_tangentspace)

    return interpolated_solution
end

function plot_interpolation_error(itp_values, true_values, parameter_points, itp_sample_parameter_point, Gr, title="")
    plot_error_interpolation = scatter(title=title, xlabel="Parameter", ylabel="Error",size=(1000,300) ,margin=5Plots.mm)
    for (index,value) in enumerate(itp_values)
        scatter!(plot_error_interpolation, [parameter_points[index]], [distance(Gr, value, true_values[index])], color = "blue", label = "", markersize = 4)
        # if parameter_points[index] in collect(itp_sample_parameter_point)
        #     scatter!(plot_error_interpolation, [parameter_points[index]], [distance(Gr, value, true_values[index])], color = "red", label = "", markersize = 6)
        # end
    end
    # add red points latter to not be obstructed b< the blue ones
    for (index,value) in enumerate(itp_values)
        # scatter!(plot_error_interpolation, [parameter_points[index]], [distance(Gr, value, true_values[index])], color = "blue", label = "", markersize = 4)
        if parameter_points[index] in collect(itp_sample_parameter_point)
            scatter!(plot_error_interpolation, [parameter_points[index]], [distance(Gr, value, true_values[index])], color = "red", label = "", markersize = 6)
        end
    end
    return plot_error_interpolation
end

function exp_log_error(points, parameter_points, origin, Gr )

    plot_error = scatter(size=(1000,300))   
    
    points_transformed = map(i -> (exp(Gr, origin  , log(Gr, origin, i ) )), points) 
    for (index,value) in enumerate(points_transformed)
        scatter!(plot_error, [[parameter_points[index]]], [distance(Gr,points[index],points_transformed[index])] , color = "blue", label = "", markersize = 5)
    end
    return plot_error
end

function plot_crossing(graphs, n_eigenfunctions, n_safety, itp_values, true_values, parameter_points, itp_sample_parameter_point, Gr)
    Σ = map(i -> hcat(Array(i.Σ)), graphs) |> x -> reduce(hcat,x) 
    plot_crossing = plot(title="$n_eigenfunctions + $n_safety Eigenvalues",xlabel="Stretch Parameter")
    for i in 1:n_eigenfunctions+n_safety
        if i <= n_eigenfunctions
            plot!(plot_crossing, parameter_points, Σ[i,:], label="",color="green")
        else
            plot!(plot_crossing, parameter_points, Σ[i,:], label="",color="red")
        end
    end
    plot_error_interpolation = plot_interpolation_error(itp_values, true_values, parameter_points, itp_sample_parameter_point, Gr, "interpolation error $n_eigenfunctions eigenfunctions")
    return plot(plot_error_interpolation, plot_crossing,layout=(2,1),size=(1000,600),margin=5Plots.mm)
end



# Custom Interpolator is used for FLUX models since the Interpolations.jl is not supported
struct CustomInterpolator
    x::Vector
    y::Vector
    
    function CustomInterpolator(y,x)
        @assert issorted(x)
        return new(x,y)
    end
end

function interpolate(citp::CustomInterpolator, x::Number)
    left_value, right_value = 0, 0
    left_index, right_index = 1, 1

    # check bound
    if x > citp.x[end] || x < citp.x[1]
        @error "Out of bounds"
        throw(DomainError(x))
    end
    
    #find the right indices 
    for (i,v) in enumerate(citp.x)
        if left_value > x
            right_value =
            right_index = i
            break
        end
        left_value = v
        left_index = i
    end

    # do a linear inter interpolation between the two selected indices
    interpolated_value = (1 - (x - left_value)/(right_value - left_value)) * citp.y[left_index] + (x - left_value)/(right_value - left_value) * citp.y[right_index]
    return interpolated_value
end
mutable struct AdaptiveInterpolator
    U::Vector{Matrix}
    parameter_space::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    interpolator::Any

    function AdaptiveInterpolator(U,parameter_space) 
        tmp  = new(U, parameter_space,  nothing)
        initialize(tmp)
        return tmp
    end
end


# #version using CustomInterpolator
# function initialize(adaptiveInterpolator::AdaptiveInterpolator)
#     U = adaptiveInterpolator.U
#     Gr = Grassmann(size(U[1])[1], size(U[1])[2])
#     parameter_space = adaptiveInterpolator.parameter_space
#     interpolator = Array{CustomInterpolator}(undef, size(parameter_space)[1])

#     for (index, value) in enumerate(U)
#         try
#             if index == 1
#                 origin = U[index]
#                 V = [log(Gr,origin,U[index]), log(Gr,origin,U[index+1])]
#                 interpolator[index] = CustomInterpolator(V,parameter_space[index:index+1])
#             end
#             if index != 1 && index != size(U)[1] 
#                 origin = U[index]
#                 V = [log(Gr,origin,U[index-1]), log(Gr,origin,U[index]), log(Gr,origin,U[index+1])]
#                 interpolator[index] = CustomInterpolator(V, parameter_space[index-1:index+1])
#             end
#             if index == size(parameter_space)[1]
#                 origin = U[index]
#                 V = [log(Gr,origin,U[index-1]), log(Gr,origin,U[index])]
#                 interpolator[index] = CustomInterpolator(V, parameter_space[index-1:index])
#             end
#         catch
#             @warn index
#         end
#     end
#     adaptiveInterpolator.interpolator = interpolator
# end

# version using Interpolations.jl
function initialize(adaptiveInterpolator::AdaptiveInterpolator)
    U = adaptiveInterpolator.U
    Gr = Grassmann(size(U[1])[1], size(U[1])[2])
    parameter_space = adaptiveInterpolator.parameter_space
    interpolator = Array{Interpolations.ScaledInterpolation}(undef, size(parameter_space)[1])

    for (index, value) in enumerate(U)
        try
            if index == 1
                origin = U[index]
                V = [log(Gr,origin,U[index]), log(Gr,origin,U[index+1])]
                interpolator[index] = interpolate(V, BSpline(Linear())) |> i->scale(i, parameter_space[index:index+1])
            end
            if index != 1 && index != size(U)[1] 
                origin = U[index]
                V = [log(Gr,origin,U[index-1]), log(Gr,origin,U[index]), log(Gr,origin,U[index+1])]
                interpolator[index] = interpolate(V, BSpline(Linear())) |> i->scale(i, parameter_space[index-1:index+1])
            end
            if index == size(parameter_space)[1]
                origin = U[index]
                V = [log(Gr,origin,U[index-1]), log(Gr,origin,U[index])]
                interpolator[index] = interpolate(V, BSpline(Linear())) |> i->scale(i, parameter_space[index-1:index])
            end
        catch
            @warn index
        end
    end
    adaptiveInterpolator.interpolator = interpolator
end

function interpolate_tangentspace(adaptiveInterpolator::AdaptiveInterpolator, parameter::Number)
    U = adaptiveInterpolator.U
    Gr = Grassmann(size(U[1])[1], size(U[1])[2])
    i = findmin(abs.(collect(adaptiveInterpolator.parameter_space).-parameter))[2]
    V = adaptiveInterpolator.interpolator[i](parameter)
    return V 

end

function interpolate(adaptiveInterpolator::AdaptiveInterpolator, parameter::Number)
    U = adaptiveInterpolator.U
    Gr = Grassmann(size(U[1])[1], size(U[1])[2])
    i = findmin(abs.(collect(adaptiveInterpolator.parameter_space).-parameter))[2]
    
    V = adaptiveInterpolator.interpolator[i](parameter)
    # V = interpolate(adaptiveInterpolator.interpolator[i],parameter)

    # compute exp function manually, because else zygote can compute derivative
    d = svd(V)
    z = U[i] * d.V * Diagonal(cos.(d.S)) * d.Vt + d.U * Diagonal(sin.(d.S)) * d.Vt
    tmp = svd(z).U
    return  tmp

    # result = exp(Gr,U[i], V)
    # return result

end

function lagrange_basis(interpolation_param, t, j)
    n = size(interpolation_param)[1]
    val = 1
    for i in 1:n
        if (i != j)
            val *= (t - interpolation_param[i] )/ ( interpolation_param[j] - interpolation_param[i] )
        end
    end 
    return val
end


