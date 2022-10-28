using LinearAlgebra
using Manifolds
using Plots
using Graphs, SimpleWeightedGraphs
using WriteVTK
using VTKDataTypes, VTKDataIO
using ManifoldsBase
using Interpolations
using PGFPlotsX
using Arpack


include("graphLaplacianToolbox.jl")
include("interpolationToolbox.jl")

#      Interpolation in normal Coordinates
#=
1. Map the dataset (eigenvectors) to the tangentspace of one of the Data points
2. Interpolate on the tangentspace, which is a R^(n-1) space
3. Map the result back to the manifold
=#


## create Data points for Interpolation
# first only a reduced order model with a stretch parameter that stretches one coordinate of the graph

system_path = "./code/Graph Laplacian/"


graphs =  CustomGraph[]
trafo_param_space = 1:10
for index in 1:size(trafo_param_space)[1]
    graphs = vcat(graphs,CustomGraph(system_path * "Meshes/spidy_$index.vtk"))
    loadMesh(graphs[index])
    createLaplacian(graphs[index])
    @assert is_connected(graphs[index].g)
    @time computeEigenGraph(graphs[index],100)
    println(index)
end 

## just need the eigenfunctions for the interpolation
n_eigenfunctions = 15
U = Array{Float64}[]
map(i -> push!(U,Matrix{Float64}(i.U[:,1:n_eigenfunctions])), graphs)


## check for crossings

n_eigenfunctions = 10  
U = Array{Float64}[] 
map(i -> push!(U,Matrix{Float64}(i.U[:,1:n_eigenfunctions])), graphs)
itp_sample_points = [U[1],U[end]]
# itp_sample_points_parameter  = 0:0.9:1
itp_sample_points_parameter  = 1:9:10
origin = itp_sample_points[1] 

interpolated_solution = interpolate_complete(itp_sample_points, itp_sample_points_parameter, trafo_param_space, origin, n_eigenfunctions)

n_safety = 1
U_safety = Array{Float64}[] 
map(i -> push!(U_safety, Matrix{Float64}(i.U[:, 1:(n_eigenfunctions+n_safety)])), graphs)
itp_sample_points_safety = [U_safety[1],U_safety[end]] 
origin_safety = itp_sample_points_safety[1] 
interpolated_solution_safety = interpolate_complete( itp_sample_points_safety, itp_sample_points_parameter, trafo_param_space, origin_safety, n_eigenfunctions + n_safety )


n = size(graphs[1].U)[1]
p = size(U[1])[2]
Gr = Grassmann(n,p)
plot_itp_error = plot_interpolation_error(interpolated_solution, U , trafo_param_space, itp_sample_points_parameter, Gr,"interpolation error $n_eigenfunctions eigenfunctions")

plot_crossing(graphs, n_eigenfunctions, n_safety,  interpolated_solution, U , trafo_param_space, itp_sample_points_parameter, Gr)


# ## 1. map data points to the tangentspace of the first data point (later it might be a better idea to choose a different "origin")
# n = size(graphs[1].U)[1]
# p = size(U[1])[2]
# Gr = Grassmann(n,p)
# interpolation_points = [U[1],U[end]]
# interpolation_param  = [1,10]
# origin = interpolation_points[1] # point where the tangentspace gets constructed 
# V = map(i -> log(Gr, origin, i), interpolation_points) # eigenvectors transformed to the tangentspace of the origin




# ## 2. interpolate on the tangentspace
# #=
# Each colum of V will be interpolated independent of the other colum
# =#
# itp = interpolate(V, BSpline(Linear()))
# sitp = scale(itp,1:9:10)    # TODO scaling function acting weird
# interpolated_tangentspace = Matrix{Float64}[]
# interpolated_tangentspace = [sitp(i) for i in trafo_param_space]
# map(i -> sitp(i), trafo_param_space)

# ## 3. retract from the tangentspace back to the manifold
# interpolated_solution = map(i -> exp(Gr,origin,i), interpolated_tangentspace)


## Compare the computed result to the actual solution
plot_error_interpolation = scatter(title="interpolation error $n_eigenfunctions eigenfunctions", xlabel="time", ylabel="Error")
for (index,value) in enumerate(interpolated_solution)
    scatter!(plot_error_interpolation, [trafo_param_space[index]], [distance(Gr, value, U[index])], color = "blue", label = "", markersize = 7)
    if index % 9 == 1
        scatter!(plot_error_interpolation, [trafo_param_space[index]], [distance(Gr, value, U[index])], color = "red", label = "", markersize = 10)
    end
end
plot_error_interpolation

## plot the eigenvalue Crossing
Σ = map(i -> hcat(Array(i.Σ)), graphs) |> x -> reduce(hcat,x) 

plot_crossings = plot(title="$n_eigenfunctions + 3 Eigenvalues",xlabel="Strech Parameter")
for i in 1:n_eigenfunctions+3
    if i <= n_eigenfunctions
        plot!(plot_crossings,trafo_param_space, Σ[i,:], label="",color="green")
    else
        plot!(plot_crossings,trafo_param_space, Σ[i,:], label="",color="red")
    end
end
plot(plot_error_interpolation,plot_crossings,layout=(2,1),size=(1000,900))


## save the Results
# for index in trafo_param_space
#     saveResults(graphs[index], system_path*"Results/"*"spidy_$index")
# end

# for index in trafo_param_space
#    graphs[index].U = interpolated_solution[index]
#     saveResults(graphs[index], system_path*"Results/"*"spidy_inter_$index")
# end

## trying to fix the eigenvalue Crossing

n_safety = 70
U_safety = Array{Float64}[]
map(i -> push!(U_safety, Matrix{Float64}(i.U[:, 1:n_eigenfunctions+n_safety])), graphs)
Gr_safety = Grassmann(n, p + n_safety)

# map to the tangentspace
interpolation_points_safety = [U_safety[1], U_safety[end]]
origin_safety = interpolation_points_safety[1] # point where the tangentspace gets constructed 
V_safety = Matrix{Float64}[]  # eigenvectors transformed to the tangentspace of the origin
map(evec -> push!(V_safety, log(Gr_safety, origin_safety, evec)), interpolation_points_safety)

# interpolate on the tangentspace
itp_safety = interpolate(V_safety, BSpline(Linear()))
sitp_safety = scale(itp_safety, 1:9:10)
interpolated_tangentspace_safety = Matrix{Float64}[]
interpolated_tangentspace_safety = [sitp_safety(i) for i in trafo_param_space]

## 3. retract from the tangentspace back to the manifold
interpolated_solution_safety = map(i -> exp(Gr_safety, origin_safety, i), interpolated_tangentspace_safety)
for (index,value) in enumerate(interpolated_solution_safety)
    scatter!(plot_error_interpolation, [trafo_param_space[index]], [distance(Gr_safety, value, U_safety[index])], color = "green", label = "", markersize = 7)
    if index % 9 == 1
        scatter!(plot_error_interpolation, [trafo_param_space[index]], [distance(Gr_safety, value, U_safety[index])], color = "red", label = "", markersize = 10)
    end
end
plot_error_interpolation

plot(plot_error_interpolation,plot_crossings,layout=(2,1),size=(1000,900))





## Try to reduce the interpolation error

# Plots.scalefontsizes(0.7)
# plots = Plots.Plot{Plots.GRBackend}[]
# for n_trial in 1:10
#     distance(Gr,interpolated_solution[n_trial],U[n_trial ])     # value to beat

#     plot_opt = scatter(xlabel="% eigenfunction $n_eigenfunctions",ylabel="% eigenfunction $(n_eigenfunctions+1)", title="Parameter t = $n_trial")
#     for i in 1:-0.1:0
#         U_opt = hcat(interpolated_solution[n_trial][:,1:end-1], i * interpolated_solution_safety[n_trial][:,n_eigenfunctions] +  (1-i) * interpolated_solution_safety[n_trial][:,n_eigenfunctions + 1])
#         scatter!(plot_opt, (i ,1-i ),marker_z= distance(Gr,U_opt,U[n_trial]), label="", c = cgrad(:thermal,rev=true))
#     end
#     push!(plots, plot_opt)
# end
# Plots.scalefontsizes()
# plot(plots[1], plots[3], plots[5], plots[6], plots[7], plots[10], size=(1000,500))




# ## test Safe Interpolator
# include("interpolationToolbox.jl")

# sft_itp = SafeInterpolator(sitp, sitp_safety)
# error_optimal = Tuple[]
# error_default= Tuple[]
# for i in 1:1:10
#     tmp = interpolate(sft_itp,i)
#     # println(size(tmp))
#     push!(error_optimal, (i, distance(Gr, tmp, U[i])))
#     default_interpolation = exp(Gr, origin,sitp(i))
#     push!(error_default, (i,distance(Gr, default_interpolation, U[i])))
# end
# plot_optimal = plot(error_optimal, seriestype=:scatter, label="crossing algorithm #1", color="blue", markersize=7, alpha = 0.5, legend=:topleft, size=(1200, 400))
# plot!(plot_optimal, error_default, seriestype=:scatter, label="default interpolation $n_eigenfunctions eigenvectors", color="red", markersize=7, alpha = 0.5)



## optimal interpolation algorithm  2 
error_rayleigh_ritz = []
for n_trial in 1:1:10
    A_tilde = interpolated_solution_safety[n_trial]' * graphs[n_trial].A * interpolated_solution_safety[n_trial]
    u_tilde = eigen(A_tilde).vectors
    u_opt = (interpolated_solution_safety[n_trial] * u_tilde)[:,1:n_eigenfunctions]

    push!(error_rayleigh_ritz, distance(Gr, U[n_trial], u_opt))
end

plot_rayleigh_ritz = plot(error_rayleigh_ritz,seriestype=:scatter, label="optimal", color="yellow", markersize=7)
plot!(plot_optimal, error_rayleigh_ritz,seriestype=:scatter, label="crossing algorithm #2 $n_safety margin eigenvectors", color="yellow", markersize=7, alpha = 0.5)
title!(plot_optimal, "Interpolation Error")
