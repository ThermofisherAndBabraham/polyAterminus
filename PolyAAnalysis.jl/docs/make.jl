using Documenter, PolyAAnalysis

makedocs(
    modules = [PolyAAnalysis],
    format = :html,
    clean = false,
    sitename = "PolyAAnalysis.jl",
    authors = "GA, KM, GM, and contributors.",
    pages = ["Methods" => "index.md"]
)

deploydocs(
    repo   = "github.com/USER/PKG.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
