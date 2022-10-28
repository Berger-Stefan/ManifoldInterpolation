using Documenter
using ManifoldInterpolation

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "ManifoldInterpolation",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [ManifoldInterpolation]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/Berger-Stefan/ManifoldInterpolation/"
)
