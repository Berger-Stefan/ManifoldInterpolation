using Documenter
using ManifoldInterpolation

makedocs(
    sitename = "ManifoldInterpolation",
    format = Documenter.HTML(),
    modules = [ManifoldInterpolation]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/Berger-Stefan/ManifoldInterpolation"
)
