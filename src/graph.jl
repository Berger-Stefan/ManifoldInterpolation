using Graphs, SimpleWeightedGraphs
using WriteVTK
using VTKDataTypes, VTKDataIO
using LinearAlgebra
using Plots
using SparseArrays
using LaTeXStrings

include("GraphLaplacianToolbox.jl")

system_path = "./code/Graph Laplacian/"

## compute the eigenfunction distance for multiple scaling

f(p,arg::Number) = diagm([1,arg,0]) * p  # scale y-direction
norm_eigenfunctions = zeros(3)
graphs = [nothing]
p = plot()
j = 2
param_space =  0.2:0.1:3


for i in param_space
    graph_tmp = CustomGraph(system_path * "Meshes/squareHighRes.vtk")
    graphs = vcat(graphs,graph_tmp)
    loadMesh(graphs[j])
    transformMesh(graphs[j],f,i)
    createLaplacian(graphs[j])
    computeEigenGraph(graphs[j])
    norm_eigenfunctions = hcat(norm_eigenfunctions, distanceEigenfunctions(graphs[2],graphs[j]))

    plot!(p,abs.([graphs[2].U[:,i]' * graphs[j].U[:,i] for i in 1:50]),label="p=$i")
    j = j + 1
end 
norm_eigenfunctions = norm_eigenfunctions[:,2:end]
graphs = graphs[2:end]


##  plotting the distance of the untransformed and transformed eigenvectors

# plot(norm_eigenfunctions[1,:], label = L"\|\| v_{p=1} - v_{p=i} \|\| ") 

saveResults(graphs[9], system_path * "Results/" * "node_08")
saveResults(graphs[11], system_path * "Results/" * "node_10")
## saving results

saveResults(graphs[1], system_path * "Results/" * "un-transformed")
saveResults(graphs[end], system_path * "Results/" * "transformed")


plot(graphs[7].Σ[1:100])
plot(graphs[9].Σ[1:100])

# ### Trafo
# # f(p) = p  # no trafo
# f(p,arg::Number) = diagm([1,arg,0]) * p  # scale y-direction
# # f(p) = map(p-> exp(p), p)
# # f(p) = map(p-> p +  0.3* rand(3), p)



# graph_untransformed = CustomGraph(system_path * "Meshes/plane_uniform.vtk")
# loadMesh(graph_untransformed)
# createLaplacian(graph_untransformed)
# computeEigenGraph(graph_untransformed)

# graph_transformed = CustomGraph(system_path * "Meshes/plane_uniform.vtk")
# loadMesh(graph_transformed)
# transformMesh(graph_transformed,f,2)
# createLaplacian(graph_transformed)
# computeEigenGraph(graph_transformed)

# ##

# saveResults(graph_transformed, system_path * "Results/" * "transformed")
# saveResults(graph_untransformed, system_path * "Results/" * "un-transformed")

# ## eigenvectors
# plot([graph_transformed.U[:,i]' * graph_untransformed.U[:,i] for i in 1:50])


# ## eigenvalues
# p = plot()
# plot!(p,graph_transformed.Σ[:],label="transformed")
# plot!(p,graph_untransformed.Σ[:],label="un-transformed")


# ## eigenfunctions 
# p = plot()
# plot!(p,graph_transformed.U[:,1],label="transformed")
# plot!(p,graph_untransformed.U[:,1],label="un-transformed")


# ##  adjencency Matrix

# heatmap(Matrix(graph_transformed.g.weights))
# heatmap(Matrix(graph_untransformed.g.weights))