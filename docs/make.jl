#!/usr/bin/env julia

using Documenter
using BenchPhaseFlow

const DOCS_DIR = @__DIR__

friendly_title(path) = begin
    base = splitext(basename(path))[1]
    title = replace(base, r"[-_]+" => " ")
    return uppercasefirst(title)
end

function collect_results_pages()
    root = joinpath(DOCS_DIR, "src", "results")
    isdir(root) || return Any[]

    by_cat = Dict{String,Vector{Pair{String,String}}}()
    for (dir, _, files) in walkdir(root)
        for file in files
            endswith(file, ".md") || continue
            abs_path = joinpath(dir, file)
            rel = relpath(abs_path, joinpath(DOCS_DIR, "src"))
            comps = splitpath(rel)
            cat = length(comps) >= 2 ? comps[2] : "results"
            entries = get!(by_cat, cat, Pair{String,String}[])
            push!(entries, friendly_title(rel) => rel)
        end
    end

    cats = sort!(collect(keys(by_cat)))
    return [cat => sort!(by_cat[cat]; by = x -> lowercase(first(x))) for cat in cats]
end

pages = Any[
    "Overview" => "index.md",
    "Articles" => [
        "BenchPhaseFlow Articles" => "articles/benchphaseflow-article1.md",
    ],

]

results_pages = collect_results_pages()
isempty(results_pages) || push!(pages, "Results" => results_pages)

html = Documenter.HTML(
    size_threshold=nothing,
    size_threshold_warn=250_000,
    inventory_version="0.1.0",
    repolink="https://github.com/PenguinxCutCell/BenchPhaseFlow.jl",
)

makedocs(
    sitename = "BenchPhaseFlow Results",
    modules = [BenchPhaseFlow],
    pages = pages,
    clean = true,
    format = html,
    repo = Documenter.Remotes.GitHub("PenguinxCutCell", "BenchPhaseFlow.jl"),
)

deploydocs(
    repo = "github.com/PenguinxCutCell/BenchPhaseFlow.jl.git",
    devbranch = "main",
    push_preview = true,
)
