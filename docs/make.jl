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
        "Quick Start" => "quickstart.md",

        "Manual" => [
            "Models" => "manual/models.md",
            "Initial States" => "manual/states.md",
            "DTWA" => "manual/dtwa.md",
            "CTWA" => "manual/ctwa.md",
            "Observables" => "manual/observables.md",
        ],

        "Examples" => [
            "Long-range Ising" => "examples/long_range_ising.md",
        ],

        "API Reference" => "api.md",
        "References & Citation" => "references.md",
    ],
)