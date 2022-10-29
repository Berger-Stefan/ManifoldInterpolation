"""

This modules allows the user to interpolate eigenvectors on the Grassmann manifold. Additionally, functions for the application of the Laplace operator on graphs is implemented. 

"""
module ManifoldInterpolation

using PyCall
using Graphs, SimpleWeightedGraphs
using VTKDataIO
using LinearAlgebra
using Plots
using SparseArrays
using LaTeXStrings
using Interpolations
using Arpack
using Manifolds
import Interpolations.interpolate

include("interpolationToolbox.jl")
include("graphLaplacianToolbox.jl")

export CustomGraph
export loadMesh
export transformMesh
export createLaplacian
export computeEigenGraph
export saveResults
export distanceEigenfunctions

export interpolate_complete
export plot_interpolation_error
export exp_log_error
export plot_crossing

end # module
