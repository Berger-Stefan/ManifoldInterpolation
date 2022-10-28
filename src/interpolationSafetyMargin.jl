# using LinearAlgebra
# using Manifolds
# using Plots
# using Graphs, SimpleWeightedGraphs
# using WriteVTK
# using VTKDataTypes, VTKDataIO
# using ManifoldsBase
# using Interpolations
# using PGFPlotsX

# include("graphLaplacianToolbox.jl")

# #=  Eigenvalue Crossing creates Problems during the interpolation SimpleWeightedGraphs

# =#
# ## create Datapoints for Interpolation

# system_path = "./code/Graph Laplacian/"

# f(p, arg::Number) = diagm([1, arg, 0]) * p  # scale y-direction
# trafo_param_space = [i for i in range(0.2, step=0.01, stop=3)]

# graphs = CustomGraph[]
# for (index, value) in enumerate(trafo_param_space)
#     graph_tmp = CustomGraph(system_path * "Meshes/square.vtk")
#     graphs = vcat(graphs, graph_tmp)
#     loadMesh(graphs[index])
#     transformMesh(graphs[index], f, value)
#     createLaplacian(graphs[index])
#     computeEigenGraph(graphs[index])
# end

# # just need the eigenfunctions for the interpolation
# n_eigenfunctions = 6
# n_safety = 2
# U = Array{Float64}[]
# U_safety = Array{Float64}[]
# map(i -> push!(U, Matrix{Float64}(i.U[:, 1:n_eigenfunctions])), graphs)
# map(i -> push!(U_safety, Matrix{Float64}(i.U[:, 1:n_eigenfunctions+n_safety])), graphs)

# ## 1. map datapoints to the tangenspace of the origin 
# n = size(graphs[1].U)[1]
# p = size(U[1])[2]

# Gr = Grassmann(n, p)
# Gr_safety = Grassmann(n, p + n_safety)

# interpolation_points = U_safety[1:20:end]
# interpolation_param = 0.2:0.2:3.0
# origin = interpolation_points[end] # point where the tangetspace gets connstructed 
# V = Matrix{Float64}[]  # eigenvectors transformed to the tangenspace of the origin
# map(evec -> push!(V, log(Gr_safety, origin, evec)), interpolation_points)

# ## 2. interpolate on the tangentspace
# #=
#     Each collum of V will be interpolated indipendent of the other collums
# =#
# itp = interpolate(V, BSpline(Linear()))
# sitp = scale(itp, interpolation_param)
# interpolated_tangentspace = Matrix{Float64}[]
# interpolated_tangentspace = [sitp(i) for i in trafo_param_space]

# ## 3. retract from the tangentspace back to the manifold
# interpolated_solution = map(i -> exp(Gr_safety, origin, i), interpolated_tangentspace)

# ## Compare the computed result to the actual solution
# # plot_error_interpolation = scatter(title="interpolation error", xlabel="Strech Parameter", ylabel="Error")
# # for (index, value) in enumerate(interpolated_soltution)
# #     scatter!(plot_error_interpolation, [trafo_param_space[index]], [distance(Gr_safty, value, U_safty[index])], color="blue", label="", markersize=4)
# #     if index % 20 == 1
# #         scatter!(plot_error_interpolation, [trafo_param_space[index]], [distance(Gr_safty, value, U_safty[index])], color="red", label="", markersize=6)
# #     end
# # end
# # plot_error_interpolation

# ## attempt to fix the interpolation between 0.4 and 0.6 (corresponds with indices 20-40)
# #=
# - inteprolat in the tangespace using the U with safty margin
# - project the graph laplacian at the interpolation point into the 

# =#

# ## How to check wheter a crossing ocours?
# # to check for a crossing, the dot products eigenfunctions of the interpolations refernce points can be computed

# abs(dot(U_safety[61][:,end-n_safety], U_safety[81][:,end-n_safety]))
# abs(dot(U_safety[81][:,end-n_safety], U_safety[101][:,end-n_safety]))

# ##
# # current problem:  when interpolating, the order of the eigenvectors can get shuffled arround 
# # question:         how to find the "right" order

# n_trial = 75 

# D_inter = exp(Gr_safety, origin, interpolated_tangentspace[n_trial])
# A_inter = graphs[n_trial].A 

# distance(Gr, D_inter[:,1:n_eigenfunctions], U[n_trial])           # distance first n eigenfunctions without safety 
# distance(Gr_safety, D_inter, U_safety[n_trial])                   # distance first n eigenfunctions with safety 

# # find D_ optimal by shuffleing eigenvector collums to reduce error
# # D_opt = hcat(D_inter[:,1],D_inter[:,3],D_inter[:,2],D_inter[:,4], D_inter[:, 5], D_inter[:,7])
# D_opt = hcat(D_inter[:,1:end], )
# distance(Gr, D_opt, U[n_trial])                     # distance with swaped eigenfunctions

# fit_inter = map(i -> abs(dot(D_opt[:,i], U[n_trial][:,i])), 1:size(D_opt)[2])



# # one problem that arrises, is that there 


# A_proj = D_inter' * Matrix{Float64}(A_inter) * D_inter
# eigenobj_inter = eigen(A_proj)
# eigenobj_inter.values
# eigenobj_inter.vectors

# eigenobj = eigen(Matrix(graphs[n_trial].A))

# eigenobj_inter.vectors




# ## save the results for comparison in paraview
# graph_inter = graphs[30]
# graph_inter.U = D_inter