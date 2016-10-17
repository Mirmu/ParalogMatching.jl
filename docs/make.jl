using Documenter, ParalogMatching

makedocs()

deploydocs(
    deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo   = "github.com/Mirmu/ParalogMatching.jl.git",
    julia  = "0.5"
)

