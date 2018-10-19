push!(LOAD_PATH, "../src")
using Documenter, PolyAAnalysis

makedocs(
    modules = [PolyAAnalysis],
    format = :html,
    clean = false,
    sitename = "PolyAAnalysis.jl",
    authors = "GA, KM, GM, and contributors.",
    pages = ["Methods" => "index.md",
             "Manual" => Any["man/guide.md",
                             "man/AnnotatePolyA.md",
                             "man/TrimmPolyA.md"
                             ]
             ]
)
