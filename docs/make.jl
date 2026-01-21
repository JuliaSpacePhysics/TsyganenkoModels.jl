using TsyganenkoModels
using Documenter

DocMeta.setdocmeta!(TsyganenkoModels, :DocTestSetup, :(using TsyganenkoModels); recursive=true)

makedocs(;
    modules=[TsyganenkoModels],
    authors="Beforerr <zzj956959688@gmail.com> and contributors",
    sitename="TsyganenkoModels.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaSpacePhysics.github.io/TsyganenkoModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSpacePhysics/TsyganenkoModels.jl",
    push_preview = true
)
