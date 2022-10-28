using LinearAlgebra
using Manifolds
using Plots
using Graphs, SimpleWeightedGraphs
using WriteVTK
using VTKDataTypes, VTKDataIO
using ManifoldsBase
using Interpolations

include("graphLaplacianToolbox.jl")
include("interpolationToolbox.jl")

function rayleigh_quotient(x, index)
    return (x' * graphs[index].A * x) /(x' * x)
end


function normalize_base!(base,p=2)
    for i in 1:size(base)[2]
        if norm(base[:,i],p) != 0 
            base[:,i] = base[:,i] / norm(base[:,i],p) 
        end
    end
end

##      Interpolation in normal Coordinates
#=
1. Map the dataset (eigenvectors) to the tangentspace of one of the data points
2. Interpolate on the tangentspace, which is a R^(n-1) space
3. Map the result back to the manifold
=#

# create data points for Interpolation
# first only a reduced order model with a stretch parameter that stretches one coordinate of the graph

system_path = "./code/Graph Laplacian/"

f(p,arg::Number) = diagm([1,arg,0]) * p  # scale y-direction
trafo_param_space = [i for i in range(0.2, step=0.01, stop=3)]

graphs =  CustomGraph[]
for (index,value) in enumerate(trafo_param_space)
    graph_tmp = CustomGraph(system_path * "Meshes/square.vtk")
    graphs = vcat(graphs,graph_tmp)
    loadMesh(graphs[index])
    transformMesh(graphs[index],f,value)
    createLaplacian(graphs[index])
    computeEigenGraph(graphs[index])
end 

##

n_eigenfunctions = 6    # crossing at 6
U = Array{Float64}[] 
map(i -> push!(U,Matrix{Float64}(i.U[:,1:n_eigenfunctions])), graphs)
itp_sample_points = U[1:20:end]
itp_sample_points_parameter  = trafo_param_space[1]:0.2:trafo_param_space[end]
origin = itp_sample_points[7] 

interpolated_solution = interpolate_complete(itp_sample_points, itp_sample_points_parameter, trafo_param_space, origin, n_eigenfunctions )

n_safety = 1
U_safety = Array{Float64}[] 
map(i -> push!(U_safety, Matrix{Float64}(i.U[:, 1:(n_eigenfunctions+n_safety)])), graphs)
itp_sample_points_safety = U_safety[1:20:end]
origin_safety = itp_sample_points_safety[end] 
interpolated_solution_safety = interpolate_complete( itp_sample_points_safety, itp_sample_points_parameter, trafo_param_space, origin_safety, n_eigenfunctions + n_safety )


Gr = Grassmann(size(U[1])[1], n_eigenfunctions)
Gr_safety = Grassmann(size(U_safety[1])[1], n_eigenfunctions + n_safety)




## Compare the computed result to the actual solution
plot_itp_error = plot_interpolation_error(interpolated_solution, U , trafo_param_space, itp_sample_points_parameter, Gr,"interpolation error $n_eigenfunctions eigenfunctions")


# # output for plots
# data = zeros(size(trafo_param_space)[1],2)
# Σ = map(i -> hcat(Array(i.Σ)), graphs) |> x -> reduce(hcat,x) 
# for (index,value) in enumerate(interpolated_solution)
#     data[index,1] = trafo_param_space[index]
#     data[index,2] = norm(value*value' - U[index]*U[index]')
#     # if trafo_param_space[index] in collect(itp_sample_points_parameter)
#         # data[index,1] = value
#         # data[index,2:end] = Σ[1:7,index]'
#     # end
#     # scatter!(plot_error_interpolation, [parameter_points[index]], [distance(Gr, value, true_values[index])], color = "blue", label = "", markersize = 4)
#     # if parameter_points[index] in collect(itp_sample_parameter_point)
#     #     scatter!(plot_error_interpolation, [parameter_points[index]], [distance(Gr, value, true_values[index])], color = "red", label = "", markersize = 6)
#     # end
# end
# using DelimitedFiles
# writedlm("crossingError.csv", data, ",")



# plot_interpolation_error([i[:,1:n_eigenfunctions] for i in interpolated_solution_safety], U , trafo_param_space, itp_sample_points_parameter, Gr_safety,"interpolation error $n_eigenfunctions eigenfunctions + $n_safety eigenfunctions")


## plot the eigenvalue Crossing
plot_crossing(graphs, n_eigenfunctions, n_safety,  interpolated_solution, U , trafo_param_space, itp_sample_points_parameter, Gr)

## crossing rayleigh ritz
# how much of Gr is in Gr_safety

# n_trial = 110
# plot_space = plot(ylabel="Similarity", xlabel="Safety Margin Size", title="Gr and Gr_safety Similarity")
# for i in n_eigenfunctions:1:size(U[1])[1] - n_eigenfunctions -50
   
#     # compute the interpolated subspaces with safety margin
#     n_safety = i
#     U_safety = Array{Float64}[] 
#     map(j -> push!(U_safety, Matrix{Float64}(j.U[:, 1:i])), graphs)
#     itp_sample_points_safety = U_safety[1:20:end]
#     origin_safety = itp_sample_points_safety[end] 
#     interpolated_solution_safety = interpolate_complete( itp_sample_points_safety, itp_sample_points_parameter, trafo_param_space, origin_safety, i )
#     # check the Similarity by doing dot dot product and summing the columns, if Similarity 1 then Gr \in Gr_safety
#     orthogonality = map( i -> norm(i,2) / n_eigenfunctions , eachcol(interpolated_solution_safety[n_trial]' * interpolated_solution[n_trial]))
#     similarity = sum(orthogonality)
#     scatter!(plot_space, [i], [similarity], label="")
# end
# plot_space

# ## error in the rayleigh ritz approach because the subspace spanned by interpolated_solution_safety does not necessary include interpolated_solution
 
n_safety = 1
U_safety = Array{Float64}[] 
map(i -> push!(U_safety, Matrix{Float64}(i.U[:, 1:(n_eigenfunctions+n_safety)])), graphs)
itp_sample_points_safety = U_safety[1:20:end]
origin_safety = itp_sample_points_safety[end] 
interpolated_solution_safety = interpolate_complete( itp_sample_points_safety, itp_sample_points_parameter, trafo_param_space, origin_safety, n_eigenfunctions + n_safety)

error_rayleigh_ritz = []
for n_trial in 1:size(trafo_param_space)[1]

    A_tilde = interpolated_solution_safety[n_trial]' * graphs[n_trial].A * interpolated_solution_safety[n_trial]
    u_tilde = eigen(A_tilde).vectors
    u_opt = (interpolated_solution_safety[n_trial] * u_tilde)[:,1:n_eigenfunctions]
    push!(error_rayleigh_ritz, [trafo_param_space[n_trial], distance(Gr, U[n_trial], u_opt)])
end

plot_rayleigh_ritz = scatter(title="interpolation error $n_eigenfunctions eigenfunctions + $n_safety safety", xlabel="Stretch Parameter", ylabel="Error",size=(1000,400))
for i in error_rayleigh_ritz
    scatter!(plot_rayleigh_ritz, [i[1]], [i[2]], label="", color="yellow")
end
plot_rayleigh_ritz


## Check if the interpolated_solution is a subset of the interpolated_solution_safety
case = 3
n_safety = 1
U_safety = Array{Float64}[] 
map(i -> push!(U_safety, Matrix{Float64}(i.U[:, 1:(n_eigenfunctions+n_safety)])), graphs)
itp_sample_points_safety = U_safety[1:20:end]
origin_safety = itp_sample_points_safety[end] 
interpolated_solution_safety = interpolate_complete( itp_sample_points_safety, itp_sample_points_parameter, trafo_param_space, origin_safety, n_eigenfunctions + n_safety)

error_rayleigh_ritz = []
for n_trial in 1:size(trafo_param_space)[1]

    if case == 1
        V_space = hcat(interpolated_solution[n_trial])                                          # 1
    elseif case == 2
        V_space = hcat(interpolated_solution_safety[n_trial])                                   # 2
    elseif case == 3
        V_space = hcat(interpolated_solution[n_trial], interpolated_solution_safety[n_trial])   # 3
    elseif case == 4 || case == 5
        # construct orthonormal space to the interpolation space
        gr_space = copy(interpolated_solution[n_trial])
        orthogonality = map( i -> norm(i,2) , eachcol(gr_safety_space' * gr_space))
        quotient_space = copy(interpolated_solution_safety[n_trial])
        for j in 1:size(gr_safety_space)[2]
            for i in 1:size(gr_space)[2]
                quotient_space[:,j] = quotient_space[:,j] - dot(quotient_space[:,j], gr_space[:,i]) * gr_space[:,i]
            end
        end
        quotient_space[abs.(quotient_space) .< eps(eltype(quotient_space))] .= zero(eltype(quotient_space)) # round all values close to 0     
        basis_quotient_space = svd(quotient_space).U
        A_tilde = basis_quotient_space' * graphs[n_trial].A * basis_quotient_space
        u_tilde = eigen(A_tilde).vectors
        u_opt = (basis_quotient_space * u_tilde)
        if case == 4
            V_space = hcat(interpolated_solution[n_trial], u_opt[:,1])                          # 4
        elseif case == 5
            V_space = hcat(interpolated_solution[n_trial], u_opt)                               # 5
        end
    end
    
    # construct a orthonormal basis in the V_space 
    V = svd(V_space).U

    # do the rayleigh ritz algorithm to get the best n_eigenfunctions that reduce the energy on the laplace operator
    A_tilde = V' * graphs[n_trial].A * V
    u_tilde = eigen(A_tilde).vectors
    u_opt = (V * u_tilde)[:,1:n_eigenfunctions]
    
    push!(error_rayleigh_ritz, [trafo_param_space[n_trial], distance(Gr, U[n_trial], u_opt)])
end

#= findings :
 - need to use a orthonormal basis to transfer A to A_tilde
 - case one ist just the default interpolation
 - case two is the interpol with a safety margin3
    - better performance then 1 almost everywhere, but patches where one is better
 - case three ist the best so far, uses  best from the combined spaces of interpolation and interpolation_points_safety
 - case four is better the one everywhere, which makes sense since the space has one extra dimension orthogonal to the interpolated space
    - but it seems like the crossing exists 
    - error can be reduced by using more columns of u_opt to construct V_space
- case 5 gives the same result as case 3 ( might be a cheaper alternative when only using n dim instead of all ?)

 - question ? Why not to use a higher dimension instead of a margin ?
    -> because if interested in 10 eigenfunction, you want the 10 th eigenfunction no 11th
=#


plot_rayleigh_ritz = scatter(title="interpolation error $n_eigenfunctions eigenfunctions + $n_safety safety", xlabel="Stretch Parameter", ylabel="Error",size=(1000,400), margin=5Plots.mm)
for i in error_rayleigh_ritz
    scatter!(plot_rayleigh_ritz, [i[1]], [i[2]], label="", color="yellow")
end
plot_rayleigh_ritz

##
# ylims!(plot_itp_error, (0,0.9))
# ylims!(plot_rayleigh_ritz, (0,0.9))
# plot_rayleigh_ritz_1 = plot_rayleigh_ritz
plot(plot_itp_error, plot_rayleigh_ritz, layout=(:,1),margin=5Plots.mm)


## comparison with polynomial interpolation and adaptiveInterpolator
# Gr = Grassmann(n,p)

interpolation_points = U[1:20:end]
interpolation_param  = 0.2:0.2:3.0
# interpolation_param  = LinRange(0.2,3.0,5) 
origin = interpolation_points[7] # point where the tangentspace gets constructed 




# interpolation_param = [0, 1/3, 2/3, 1]
plot_lagrange = plot(title="lagrange basis, $(size(interpolation_param)[1]) order")
for (index,value) in enumerate(interpolation_param)
    plot!(plot_lagrange, [ interpolation_param[1]:0.01:interpolation_param[end]],[lagrange_basis(interpolation_param, i, index) for i in interpolation_param[1]:0.01:interpolation_param[end]], label="$index")
end
plot_lagrange


Gr = Grassmann(size(U[1])[1], n_eigenfunctions)
lagrange_interpolation = []
for i in trafo_param_space
    interpolated_value = zeros(size(itp_sample_points[1])[1], n_eigenfunctions)
    for j in 1:size(interpolation_param)[1]
        interpolated_value += lagrange_basis(interpolation_param, i, j) .* log(Gr, origin, interpolation_points[j])
    end
    push!(lagrange_interpolation, exp(Gr, origin, interpolated_value))
end



adaptive_interpolator = AdaptiveInterpolator(interpolation_points, interpolation_param)
adaptive_interpolation = []
for i in trafo_param_space
    push!(adaptive_interpolation, interpolate(adaptive_interpolator,i))
end
plot_itp_error_adaptive = plot_interpolation_error(adaptive_interpolation, U , trafo_param_space, itp_sample_points_parameter, Gr,"interpolation error adaptive piecewise linear functions")

plot_itp_error = plot_interpolation_error(interpolated_solution, U , trafo_param_space, itp_sample_points_parameter, Gr,"interpolation error $n_eigenfunctions eigenfunctions Splines at one origin")

plot_itp_error_polynomial = plot_interpolation_error(lagrange_interpolation, U , trafo_param_space, itp_sample_points_parameter, Gr,"interpolation error lagrange polynomial")
plot_itp_error_polynomial_ylim =  ylims!(plot_itp_error_polynomial,(0,0.001))
title!(plot_itp_error_polynomial_ylim, "interpolation error lagrange polynomial (Y axis zoomed)")
plot_itp_error_polynomial = plot_interpolation_error(lagrange_interpolation, U , trafo_param_space, itp_sample_points_parameter, Gr,"interpolation error lagrange polynomial")

plot(plot_itp_error_polynomial, plot_itp_error_polynomial_ylim, plot_itp_error_adaptive, plot_itp_error , layout=(:,1), size=(1000,900))