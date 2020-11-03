using WCSRO
using Documenter

makedocs(;
    modules=[WCSRO],
    authors="Ulrich Pototschnig <ulrich.pototschnig@hotmail.com> and contributors",
    repo="https://github.com/ulipoto/WCSRO.jl/blob/{commit}{path}#L{line}",
    sitename="WCSRO.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ulipoto.github.io/WCSRO.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ulipoto/WCSRO.jl",
)
