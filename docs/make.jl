using Phase4s
using Documenter

DocMeta.setdocmeta!(Phase4s, :DocTestSetup, :(using Phase4s); recursive=true)

makedocs(;
    modules=[Phase4s],
    authors="John Lapeyre",
    sitename="Phase4s.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
