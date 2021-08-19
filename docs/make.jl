##
using Documenter,JSWAP

makedocs(
    modules=[JSWAP],
    authors="Yi Zhang",
    sitename="JuliaSoundWAvePropagation",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
        "Inputs" =>"Inputs.md"
        "Outputs"=>"Outputs.md"
        "Forward simulation"=>"Forward_simulation.md"
        "Utilities"=>"Utilities.md"
        ],
        "List of functions" =>"list_of_functions.md"
    ],
)

deploydocs(;
    repo="github.com/deconvolution/JSWAP.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    devurl = "dev",
    forcepush=true,
    push_preview =true
)
