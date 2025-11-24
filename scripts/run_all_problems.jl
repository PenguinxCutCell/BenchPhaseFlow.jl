#!/usr/bin/env julia

"""
Run every benchmark problem script in `problems/` and fail if any of them exit non-zero.

Usage
-----
    julia --project scripts/run_all_problems.jl [filters...]

- Pass one or more substring filters to run only the matching scripts.
- Add `--list` to print the current manifest and exit.
"""

using Printf
using Logging

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const JULIA_CMD = Base.julia_cmd()
const PROBLEMS_DIR = joinpath(REPO_ROOT, "problems")

function discover_problem_scripts()
    isdir(PROBLEMS_DIR) || error("Problems directory not found at $(PROBLEMS_DIR)")
    scripts = String[]
    for (root, _, files) in walkdir(PROBLEMS_DIR)
        for file in files
            endswith(file, ".jl") || continue
            abs_path = joinpath(root, file)
            push!(scripts, relpath(abs_path, REPO_ROOT))
        end
    end
    sort!(scripts)
    return scripts
end

struct ProblemResult
    path::String
    success::Bool
    elapsed::Float64
    error::Union{Nothing,Exception}
end

function format_seconds(seconds::Float64)
    return seconds < 120 ? @sprintf("%.1fs", seconds) : @sprintf("%.1fm", seconds / 60)
end

function select_scripts(args::Vector{String})
    all_scripts = discover_problem_scripts()
    filters = String[]
    list_only = false
    for arg in args
        if arg in ("--help", "-h")
            println("Usage: julia --project scripts/run_all_problems.jl [--list] [filters...]")
            println("Run all benchmark scripts, optionally filtering by substring.")
            exit(0)
        elseif arg == "--list"
            list_only = true
        else
            push!(filters, lowercase(arg))
        end
    end

    if list_only
        println("Registered problem scripts:")
        for script in all_scripts
            println("  " * script)
        end
        exit(0)
    end

    return isempty(filters) ? all_scripts : [script for script in all_scripts if any(f -> occursin(f, lowercase(script)), filters)]
end

function run_problem(script::String)
    abs_path = normpath(joinpath(REPO_ROOT, script))
    cmd = `$(JULIA_CMD) --project=$(REPO_ROOT) --color=yes --startup-file=no $(abs_path)`
    println("\n>>> Running $(script)")
    start_time = time()
    try
        run(cmd)
        elapsed = time() - start_time
        println("<<< Completed $(script) in $(format_seconds(elapsed))")
        return ProblemResult(script, true, elapsed, nothing)
    catch err
        elapsed = time() - start_time
        if err isa Base.ProcessFailedException
            @error "Script exited with non-zero status" script=script exit_code=err.cmd.process.exitcode
        else
            @error "Script run threw" script=script exception=err
        end
        return ProblemResult(script, false, elapsed, err)
    end
end

function main()
    scripts = select_scripts(copy(ARGS))
    isempty(scripts) && error("No scripts matched the provided filters")
    println("Running $(length(scripts)) problem script(s)...")

    results = ProblemResult[]
    for script in scripts
        push!(results, run_problem(script))
    end

    failures = filter(r -> !r.success, results)
    total_time = sum(r.elapsed for r in results)
    println("\n=== Problem summary ===")
    println("Total: $(length(results)) | Passed: $(length(results) - length(failures)) | Failed: $(length(failures)) | Wall time: $(format_seconds(total_time))")
    if !isempty(failures)
        println("Failed scripts:")
        for failure in failures
            println("  $(failure.path) (in $(format_seconds(failure.elapsed)))")
        end
        error("Problem scripts failed")
    end
end

main()
