using Documenter, PolyAAnalysis

makedocs(
    modules = [PolyAAnalysis],
    format = :html,
    clean = false,
    sitename = "PolyAAnalysis.jl",
    authors = "GA, KM, GM, and contributors.",
    pages = ["Methods" => "index.md",
             "Manual" => "man/guide.md"
                        
             ]
)
