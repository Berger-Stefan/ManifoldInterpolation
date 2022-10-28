mutable struct CustomGraph
	# attributes
	pathMesh::String
	N::Integer
	points 
	mesh 
	Σ::Array{Float64} 
	U 
	g 
	A 
	# default constructor
	function CustomGraph(pathMesh::String)
		return new(pathMesh,50,nothing,nothing,zeros(1),nothing,nothing,nothing)
	end
	# N number of eigenfunctions
	function CustomGraph(pathMesh::String, N::Integer)
		return new(pathMesh,N,nothing,nothing,zeros(1),nothing,nothing,nothing)
	end
end

function loadMesh(graph::CustomGraph)
	graph.mesh = read_vtk(graph.pathMesh)
	graph.points = [graph.mesh.point_coords[:,i] for i in 1:size(graph.mesh.point_coords,2)]
end

function transformMesh(graph::CustomGraph, transformationOperator::Function,functionArg)
	graph.points = transformationOperator.(graph.points,functionArg ...)
	graph.mesh.point_coords = hcat(graph.points...)	
end

function createLaplacian(graph::CustomGraph)
	# Points
	graph.points = [graph.mesh.point_coords[:,i] for i in 1:size(graph.mesh.point_coords,2)]	# write point coord from mesh to graph
	nv = length(graph.points)

	# Graph
	graph.g = SimpleWeightedGraph(nv)
	elements = filter(x->length(x)==3, graph.mesh.cell_connectivity) # filter for elements of size 3
	ne = length(elements)
	for el in 1:ne
	# Add edges
	edges = Iterators.product(elements[el], elements[el]) |> collect |> vec
	map(p -> add_edge!(graph.g, p[1], p[2], norm(graph.points[p[1]] - graph.points[p[2]])), edges)
	end

	# Graph Laplacian
	graph.A = laplacian_matrix(graph.g)
end

# function computeEigenGraph(graph::CustomGraph,N::Int )
# 	# Diagonalization
# 	graph.N = N
# 	eigobj = eigen(Array(graph.A))
# 	graph.U = eigobj.vectors[:,1:N]
# 	graph.Σ = eigobj.values[1:N]
# end

function computeEigenGraph(graph::CustomGraph,N::Int )
	# Diagonalization
	graph.N = N
	eigobj = eigs(graph.A,nev=N, which=:SM, maxiter = 1000)
	graph.U = eigobj[2]
	graph.Σ = eigobj[1]
end


function computeEigenGraph(graph::CustomGraph)
	# Diagonalization
	eigobj = eigen(Array(graph.A))
	graph.U = eigobj.vectors
	graph.Σ = eigobj.values
end


function saveResults(graph::CustomGraph, outputPath::String="result")
	# Write Output VTK
	for i=1:size(graph.U)[2]
		meshtmp = graph.mesh
		meshtmp.point_data["eigvec"] = graph.U[:,i]
		meshtmp.point_data["eigval"] = graph.Σ[i] * ones(length(graph.points))
		write_vtk(graph.mesh, outputPath * "_$i.vtu")
	end
end

function distanceEigenfunctions(graph1::CustomGraph, graph2::CustomGraph, n::Number=50)
	norm_inf =  norm(graph1.U[:,n] - graph2.U[:,n], 1)
	norm_fro = norm(graph1.U[:,n] - graph2.U[:,n], 2)
	norm_max = norm(graph1.U[:,n] - graph2.U[:,n], Inf)
	return [norm_inf, norm_fro, norm_max]
end