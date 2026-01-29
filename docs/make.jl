#!/usr/bin/env julia

using Documenter
using CSV
using DataFrames
using Printf

const DOCS_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(DOCS_DIR, ".."))
const RESULTS_DIR = joinpath(PROJECT_ROOT, "results")
const PROBLEMS_DIR = joinpath(PROJECT_ROOT, "problems")
const GENERATED_RESULTS_DIR = joinpath(DOCS_DIR, "src", "results")

mkpath(GENERATED_RESULTS_DIR)

function clean_generated_pages()
    for entry in readdir(GENERATED_RESULTS_DIR; join=true)
        endswith(entry, ".md") && rm(entry; force=true)
    end
end

function extract_docstring(text::String)
    m = match(r"\"\"\"(?s)(.*?)\"\"\"", text)
    return m === nothing ? nothing : strip(m.captures[1])
end

function collect_problem_descriptions()
    descriptions = Dict{String,String}()
    isdir(PROBLEMS_DIR) || return descriptions
    for (root, _, files) in walkdir(PROBLEMS_DIR)
        for file in files
            endswith(file, ".jl") || continue
            abs_path = joinpath(root, file)
            contents = read(abs_path, String)
            desc = something(extract_docstring(contents), "Description missing in $(relpath(abs_path, PROJECT_ROOT))")
            base = replace(file, ".jl" => "")
            descriptions[base] = desc
            for m in eachmatch(r"write_convergence_csv\(\s*\"([^\"]+)\"", contents)
                descriptions[m.captures[1]] = desc
            end
        end
    end
    return descriptions
end

function describe_result(key::String, descriptions::Dict{String,String})
    get(descriptions, key) do
        stripped = replace(key, r"_(Convergence|Diagnostics|Temporal|Mass|MeanError|Linf)$" => "")
        get(descriptions, stripped, "No description available; please add a docstring to the originating problem script.")
    end
end

function slugify(path::String)
    slug = lowercase(replace(path, r"[^A-Za-z0-9]+" => "-"))
    slug = strip(slug, ['-'])
    isempty(slug) && (slug = "result")
    return slug
end

function format_cell(value)
    if value isa AbstractFloat
        if isnan(value)
            return "NaN"
        elseif isfinite(value)
            return @sprintf("%.6g", value)
        else
            return string(value)
        end
    elseif value isa Integer
        return string(value)
    elseif value isa Bool
        return value ? "true" : "false"
    elseif value isa Missing
        return "missing"
    else
        return string(value)
    end
end

function dataframe_to_markdown(df::DataFrame)
    if nrow(df) == 0
        return "_CSV has no rows._"
    end
    header = join(string.(names(df)), " | ")
    divider = join(fill("---", ncol(df)), " | ")
    rows = String[
        join([format_cell(row[col]) for col in names(df)], " | ") for row in eachrow(df)
    ]
    table_lines = Vector{String}()
    push!(table_lines, "| " * header * " |")
    push!(table_lines, "| " * divider * " |")
    for row in rows
        push!(table_lines, "| " * row * " |")
    end
    return join(table_lines, "\n")
end

struct ResultEntry
    title::String
    nav_title::String
    csv_abs::String
    csv_rel::String
    page_file::String
    description::String
end

function collect_results(descriptions::Dict{String,String})
    entries = ResultEntry[]
    isdir(RESULTS_DIR) || return entries
    csv_paths = String[]
    for (root, _, files) in walkdir(RESULTS_DIR)
        for file in files
            endswith(file, ".csv") || continue
            push!(csv_paths, joinpath(root, file))
        end
    end
    sort!(csv_paths)
    for csv_path in csv_paths
        rel_csv = relpath(csv_path, PROJECT_ROOT)
        rel_without_ext = replace(rel_csv, r"\.csv$" => "")
        base = splitext(basename(csv_path))[1]
        description = describe_result(base, descriptions)
        title = replace(base, "_" => " ")
        relative_title = replace(rel_without_ext, "results/" => "")
        nav_title = join(split(relative_title, '/'), " / ")
        slug = slugify(rel_without_ext)
        page_file = slug * ".md"
        push!(entries, ResultEntry(title, nav_title, csv_path, rel_csv, page_file, description))
    end
    return entries
end

function write_result_page(entry::ResultEntry)
    df = CSV.read(entry.csv_abs, DataFrame)
    table_md = dataframe_to_markdown(df)
    page_path = joinpath(GENERATED_RESULTS_DIR, entry.page_file)
    open(page_path, "w") do io
        println(io, "# ", entry.title)
        println(io)
        println(io, entry.description)
        println(io)
        println(io, "**CSV source:** `", entry.csv_rel, "`")
        println(io)
        println(io, table_md)
    end
end

clean_generated_pages()
descriptions = collect_problem_descriptions()
result_entries = collect_results(descriptions)
for entry in result_entries
    write_result_page(entry)
end

pages = Any[
    "Overview" => "index.md"
]

if !isempty(result_entries)
    push!(pages, "Results" => [entry.nav_title => joinpath("results", entry.page_file) for entry in result_entries])
end

html = Documenter.HTML(
    size_threshold=nothing, 
    size_threshold_warn=250_000, 
    inventory_version="0.1.0",
    repolink="https://github.com/PenguinxCutCell/BenchPhaseFlow.jl",
)

makedocs(
    sitename = "BenchPhaseFlow Results",
    modules = Module[],
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
