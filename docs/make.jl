using Documenter
using ManifoldInterpolation

PAGES = [
            "Home" => "index.md",
            "Documentation" => ["Interpolation.md", "GraphLaplacian.md"],
        ]


makedocs(
    sitename = "ManifoldInterpolation",
    format = Documenter.HTML(prettyurls = false ),
    modules = [ManifoldInterpolation],
    # pages = PAGES
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    devbranch = "main",
    repo = "https://github.com/Berger-Stefan/ManifoldInterpolation"
)
