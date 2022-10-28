# using Graphs, SimpleWeightedGraphs, GraphPlot, GraphRecipes
# using Zygote
# using LinearAlgebra
# using DataFrames
# using Dates
# using Cairo
# using Manifolds
# using Interpolations
# using Plots
# using Colors, ColorSchemes
# using Compose

# include("interpolationToolbox.jl")

# function frobenius_norm(U_approx, U_exact)
#     norm(U_approx * U_approx' - U_exact*U_exact',2)
# end

# function main_component(g)
#     c = connected_components(g)
#     _, i = findmax(length.(c))
#     return (c[i], g[c[i]])
# end

# function add_padding(val_1, val_2)
#     eigenvec_padded = eigenvec[val_1]
#     padding = []
#     for (index,comp) in enumerate(connected_components(graphs[val_1]))
#         if index > 1
#             for i in comp
#                 if i in connected_components(graphs[val_2])[1]
#                     append!(padding, i)
#                 end
#             end
#         end
#     end
#     padding_1_sorted = sort(padding)

#     for i in padding_1_sorted
#         eigenvec_padded = vcat(eigenvec_padded[1:i-1,:], zeros(1,100), eigenvec_padded[i:end,:])
#     end

#     eigenvec_padded,padding
# end

# function remove_padding(eigenvec_padded, padding)
#     eigenvec = eigenvec_padded
#     # for i in padding
#     #     eigenvec = vcat(eigenvec[1:i-1,:], eigenvec[i+1:end,:])
#     # end
#     eigenvec[setdiff(1:end,padding),:]
# end

# function colorscale(value)
#     c = colormap("RdBu") 
#     v = floor(value * 50) 
#     return c[Int(v+50)]
# end
# ##
# # t=1
# graph_1 = SimpleWeightedGraph(7)
# add_edge!(graph_1, 1, 2)
# add_edge!(graph_1, 2, 3)
# add_edge!(graph_1, 3, 4)

# # t=3
# graph_2 = SimpleWeightedGraph(7)
# add_edge!(graph_2, 1, 2)
# add_edge!(graph_2, 2, 3)
# add_edge!(graph_2, 3, 4)
# add_edge!(graph_2, 4, 5)

# # t=3
# graph_3 = SimpleWeightedGraph(7)
# add_edge!(graph_3, 1, 2)
# add_edge!(graph_3, 2, 3)
# add_edge!(graph_3, 3, 4)
# add_edge!(graph_3, 4, 5)
# add_edge!(graph_3, 5, 6)

# ## 

# indices_1, components_1 = main_component(graph_1)
# indices_2, components_2 = main_component(graph_2)
# indices_3, components_3 = main_component(graph_3)

# eigenvectors_1 = eigen(Array(laplacian_matrix(components_1))).vectors
# eigenvectors_3 = eigen(Array(laplacian_matrix(components_3))).vectors
# eigenvectors_2 = eigen(Array(laplacian_matrix(components_2))).vectors


# gplot(graph_1, nodelabel = 1:nv(graph_1))
# gplot(graph_2, nodelabel = 1:nv(graph_2))
# gplot(graph_3, nodelabel = 1:nv(graph_3))


# ##

# n_eigenfunctions = 4 # TODO More eigenfunctions

# # find p(size biggest common component) 
# p = size(components_1)[1] > size(components_3)[1] ? size(components_1)[1] : size(components_3)[1]

# # find the common nodes
# eigenvectors_1_modified = vcat(eigenvectors_1[1:2,1:n_eigenfunctions], eigenvectors_1[3:end,1:n_eigenfunctions], zeros(2,n_eigenfunctions))
# eigenvectors_2_modified = vcat(eigenvectors_2[:,1:n_eigenfunctions], zeros(1,n_eigenfunctions))
# eigenvectors_3_modified = eigenvectors_3[:,1:n_eigenfunctions]

# Gr = Grassmann(p, n_eigenfunctions)
# origin = eigenvectors_3_modified
# interpolation_points = [eigenvectors_1_modified, eigenvectors_3_modified]
# V = map(i -> log(Gr, origin, i), interpolation_points)
# sitp = interpolate(V, BSpline(Linear())) |> i -> Interpolations.scale(i, 1:2:3)
# itp = sitp

# # TODO after adding extra rows, how to remove the rows again
# eigenvectors_2_interpolated = sitp(2) |> i-> exp(Gr, origin, i)


# Gr_reduced = Grassmann(p-1, n_eigenfunctions)
# eigenvectors_2_tangentspace = sitp(2)
# eigenvectors_2_tangentspace = vcat(eigenvectors_2_tangentspace[1:5,:], zeros(4)')
# eigenvectors_2_interpolated = exp(Gr, origin, eigenvectors_2_tangentspace)

# # use rayleigh ritz to improve the result
# V = eigenvectors_2_interpolated[1:5,:]
# V = svd(V).U


# A_tilde = V' * Array(laplacian_matrix(components_2)) * V
# u_tilde = eigen(A_tilde).vectors
# u_opt = (V * u_tilde)
# u_opt = u_opt[1:5,:]
# u_opt = svd(u_opt).U
# # TODO u_opt und eigenvectors_2_modified not comparable
# Gr_reduced = Grassmann(p-1, n_eigenfunctions)
# distance(Gr_reduced, u_opt, eigenvectors_2[:,1:n_eigenfunctions])

# ## visualization
# n_trial = 2

# values = eigenvectors_2[:,n_trial]
# nodefillc = [colorscale(i) for i in values];
# plot_cal = gplot(components_2, nodelabel = 1:nv(components_2), nodefillc = nodefillc)

# values = eigenvectors_2_interpolated[1:5,n_trial]
# nodefillc = [colorscale(i) for i in values];
# plot_int = gplot(components_2, nodelabel = 1:nv(components_2), nodefillc = nodefillc)


# values = u_opt[:,n_trial]
# nodefillc = [colorscale(i) for i in values];
# plot_opt = gplot(components_2, nodelabel = 1:nv(components_2), nodefillc = nodefillc)

# values = eigenvectors_3[:,n_trial]
# nodefillc = [colorscale(i) for i in values];
# plot_graph_3 = gplot(components_3, nodelabel = 1:nv(components_3), nodefillc = nodefillc)


# set_default_graphic_size(40cm, 20cm)
# Compose.compose( context(),
#                 (context(0.1, 0.1, 0.3, 0.05),
#                     (context(), Compose.text(0.3,0,"real Solution"), fontsize(18pt),fill("white")),
#                     (context(0,0.6),plot_cal)),
#                 (context(0.1, 0.3, 0.3, 0.05),
#                     (context(), Compose.text(0.3,0,"interpolated Solution"), fontsize(18pt),fill("white")),
#                     (context(0,0.6),plot_int)),
#                 (context(0.1, 0.5, 0.3, 0.05),
#                     (context(), Compose.text(0.3,0,"optimal Solution"), fontsize(18pt),fill("white")),
#                     (context(0,0.6),plot_opt)),
#                 (context(0.1, 0.7, 0.3, 0.05),
#                     (context(), Compose.text(0.3,0,"t=3 Solution"), fontsize(18pt),fill("white")),
#                     (context(0,0.6),plot_graph_3)))



# ##
# A = Array(laplacian_matrix(main_component(graph_3)[2]))
# # A = -1 .* A
# A[4,4] = 1
# A[5,5] += 10000
# A[6,6] += 10000
# tmp = eigen(A).vectors[:,1:n_eigenfunctions]
# frobenius_norm(tmp[1:4,:],eigenvectors_1[:,1:n_eigenfunctions])


# ## instead of cutting values attempt with potential
# gplot(graph_1, nodelabel = 1:nv(graph_1))
# gplot(graph_3, nodelabel = 1:nv(graph_3))


# potential_values = [2^i for i in 0:1:20]
# error_ = Matrix(undef, 2, size(potential_values)[1])
# value_ = Vector(undef, size(potential_values)[1])
# for (index,value) in enumerate(potential_values)
#     A = Array(laplacian_matrix(main_component(graph_3)[2]))
#     A[4,4] = 1
#     A[5,5] += value
#     A[6,6] += value 

#     n_eigenfunctions = 3
#     value_[index] = eigen(A).vectors[:,1:n_eigenfunctions]

#     tmp = value_[index][1:4,1:n_eigenfunctions]
#     er_1 = frobenius_norm(tmp,eigenvectors_1[:,1:n_eigenfunctions])
#     er_2 = distance(Grassmann(4,4),tmp, eigenvectors_1[:,1:n_eigenfunctions])

#     error_[1,index] = er_1
#     error_[2,index] = er_2
# end

# plot( potential_values, error_[1,:],marker=:auto, label ="Error Frobenius", xaxis=:log, xlabel="Potetial Size", ylabel="Error")
# plot!(potential_values, error_[2,:],marker=:auto, label ="Error Subspace Angle")


# ## find a better way then interpolation
# n_eigenfunctions  = 3
# Gr = Grassmann(6,n_eigenfunctions)

# x_start = exp(Gr,origin,itp(2))[:,1:n_eigenfunctions]
# x = x_start   

# μ = [1,0.05]
# null_rows = [6]
# function F(x::Matrix)
#     tmp1 = sum([norm(x[i,:])^2 for i in null_rows])
#     # tmp2 = (distance(Gr,x, x_start))^2
#     return μ[1] *tmp1 + μ[2] *frobenius_norm(x,x_start)
# end

# function residuum(x)
#     return abs(distance(Gr, x[1:5,1:n_eigenfunctions], eigenvectors_2[:,1:n_eigenfunctions]))
# end 

# step_size = 0.002
# # while norm(Manifolds._gradient(F, x, Manifolds.ZygoteDiffBackend()),2 ) > 0.0001
# #     println(residuum(x),"\t", 
# #             norm(Manifolds._gradient(F, x, Manifolds.ZygoteDiffBackend()),2 ) , "\t", 
# #             sum([norm(x[i,:])^2 for i in null_rows]), "\t",
# #             F(x))
# #     x = exp(Gr, x, - step_size .*  project(Gr, x, Manifolds._gradient(F, x, Manifolds.ZygoteDiffBackend())))
# # end

# result = x[1:5,1:n_eigenfunctions]

# distance(Grassmann(5,n_eigenfunctions), result, eigenvectors_2[:,1:n_eigenfunctions])
# distance(Grassmann(5,n_eigenfunctions), exp(Gr,origin,itp(2))[1:5,1:n_eigenfunctions], eigenvectors_2[:,1:n_eigenfunctions])

# ## instead potential operator V(t), add a diffusion operator 

# A = Array(laplacian_matrix(main_component(graph_3)[2]))

# # Algorithm for finding the the right "Potential" Operator 
# # - identify which nodes are not present
# # - identify to which nodes those connect and reduce the rank of the nodes by one
# # - set value of the nodes that are not present to a high value and enforce a "Potential"
# # - set the distance to the inverse potential


# A[5,5] = 1
# A[6,6] = 1e10
# A[6,5] = 1e-10
# A[5,6] = 1e-10

# n_eigenfunctions = 2 

# tmp = eigen(A).vectors[:,1:n_eigenfunctions]
# frobenius_norm(tmp[1:5,:],eigenvectors_2[:,1:n_eigenfunctions])
# # distance(Grassmann(6,3), tmp[1:5,:], eigenvectors_2[:,1:n_eigenfunctions])