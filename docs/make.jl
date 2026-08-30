using Documenter
using TWA

makedocs(
    sitename = "TWA.jl",
    modules = [TWA],
    checkdocs = :none,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
    ],
)