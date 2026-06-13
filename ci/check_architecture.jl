#!/usr/bin/env julia

using FuncDeps
using TOML

struct Violation
    rule::String
    src::String
    dst::String
    detail::String
end

function usage()
    println("""
Usage:
  julia --project=ci/FuncDepsEnv ci/check_architecture.jl <source_dir> <policy_file> <report_file>

Example:
  julia --project=ci/FuncDepsEnv ci/check_architecture.jl src ci/architecture_policy.toml artifacts/funcdeps/architecture_report.txt
""")
end

length(ARGS) >= 3 || (usage(); exit(2))

source_dir = ARGS[1]
policy_file = ARGS[2]
report_file = ARGS[3]

isdir(source_dir) || error("source_dir does not exist: $source_dir")
isfile(policy_file) || error("policy_file does not exist: $policy_file")
mkpath(dirname(report_file))

policy = TOML.parsefile(policy_file)
case_insensitive = get(policy, "case_insensitive", true)

normname(s) = begin
    t = String(s)
    t = replace(t, "\\" => "/")
    t = basename(t)
    t = replace(t, r"\.jl$" => "")
    case_insensitive ? lowercase(t) : t
end

normfull(s) = begin
    t = String(s)
    t = replace(t, "\\" => "/")
    t = replace(t, r"\.jl" => "")
    case_insensitive ? lowercase(t) : t
end

function module_norm(full_name::AbstractString)
    return normname(module_of(full_name))
end

function file_norm(path::AbstractString)
    return normname(path)
end

function string_list(x)
    x === nothing && return String[]
    x isa AbstractVector || return String[]
    return String.(x)
end

function get_info_file(idx, full_name::AbstractString)
    info = get(idx.infos, String(full_name), nothing)
    info === nothing && return ""
    return info.file
end

function edge_key(src, dst)
    string(normfull(src), " -> ", normfull(dst))
end

ignored_edges = Set{String}()
for e in string_list(get(policy, "ignored_edges", String[]))
    push!(ignored_edges, normfull(e))
end

function is_ignored(src, dst)
    k = edge_key(src, dst)
    return k in ignored_edges
end

infos = scan_project(source_dir)
idx = build_index(infos)

violations = Violation[]

# 1. independent_modules: functions in these modules may not call project
# functions from any other discovered module.
for mod in string_list(get(policy, "independent_modules", String[]))
    m = normname(mod)
    for info in infos
        module_norm(info.full_name) == m || continue
        for callee in info.calls
            module_norm(callee) == m && continue
            is_ignored(info.full_name, callee) && continue
            push!(violations, Violation(
                "independent_modules",
                info.full_name,
                callee,
                "module $(module_of(info.full_name)) is marked independent but calls module $(module_of(callee))",
            ))
        end
    end
end

# 2. independent_files: functions in these files may not call project functions
# defined in a different file.
for file in string_list(get(policy, "independent_files", String[]))
    f = file_norm(file)
    for info in infos
        file_norm(info.file) == f || continue
        for callee in info.calls
            callee_file = get_info_file(idx, callee)
            isempty(callee_file) && continue
            file_norm(callee_file) == f && continue
            is_ignored(info.full_name, callee) && continue
            push!(violations, Violation(
                "independent_files",
                info.full_name,
                callee,
                "file $(basename(info.file)) is marked independent but calls function in $(basename(callee_file))",
            ))
        end
    end
end

# 3. forbidden_module_deps: caller-module => blocked callee modules.
for (srcmod, blocked) in get(policy, "forbidden_module_deps", Dict{String,Any}())
    src_m = normname(srcmod)
    blocked_set = Set(normname.(string_list(blocked)))
    for info in infos
        module_norm(info.full_name) == src_m || continue
        for callee in info.calls
            dst_m = module_norm(callee)
            dst_m in blocked_set || continue
            is_ignored(info.full_name, callee) && continue
            push!(violations, Violation(
                "forbidden_module_deps",
                info.full_name,
                callee,
                "module $(module_of(info.full_name)) is forbidden from depending on module $(module_of(callee))",
            ))
        end
    end
end

# 4. optional allowed_module_deps: caller-module => only these cross-module
# dependencies are allowed.
for (srcmod, allowed) in get(policy, "allowed_module_deps", Dict{String,Any}())
    src_m = normname(srcmod)
    allowed_set = Set(normname.(string_list(allowed)))
    for info in infos
        module_norm(info.full_name) == src_m || continue
        for callee in info.calls
            dst_m = module_norm(callee)
            dst_m == src_m && continue
            dst_m in allowed_set && continue
            is_ignored(info.full_name, callee) && continue
            push!(violations, Violation(
                "allowed_module_deps",
                info.full_name,
                callee,
                "module $(module_of(info.full_name)) may only depend on $(sort(collect(allowed_set))) but calls $(module_of(callee))",
            ))
        end
    end
end

# 5. independent_functions.
ind_funcs = get(policy, "independent_functions", Dict{String,Any}())
if ind_funcs isa Dict
    for fname in string_list(get(ind_funcs, "cross_module", String[]))
        target = normfull(fname)
        for info in infos
            normfull(info.full_name) == target || continue
            for callee in info.calls
                module_norm(callee) == module_norm(info.full_name) && continue
                is_ignored(info.full_name, callee) && continue
                push!(violations, Violation(
                    "independent_functions.cross_module",
                    info.full_name,
                    callee,
                    "function $(info.full_name) is marked cross-module independent",
                ))
            end
        end
    end

    for fname in string_list(get(ind_funcs, "no_internal_calls", String[]))
        target = normfull(fname)
        for info in infos
            normfull(info.full_name) == target || continue
            for callee in info.calls
                is_ignored(info.full_name, callee) && continue
                push!(violations, Violation(
                    "independent_functions.no_internal_calls",
                    info.full_name,
                    callee,
                    "function $(info.full_name) is marked as having no internal project calls",
                ))
            end
        end
    end
end

# 6. exact forbidden function deps.
for (srcfun, blocked) in get(policy, "forbidden_function_deps", Dict{String,Any}())
    src_f = normfull(srcfun)
    blocked_set = Set(normfull.(string_list(blocked)))
    for info in infos
        normfull(info.full_name) == src_f || continue
        for callee in info.calls
            normfull(callee) in blocked_set || continue
            is_ignored(info.full_name, callee) && continue
            push!(violations, Violation(
                "forbidden_function_deps",
                info.full_name,
                callee,
                "exact function dependency is forbidden by policy",
            ))
        end
    end
end

open(report_file, "w") do io
    println(io, "FusionRings.jl architecture report")
    println(io, "===================================")
    println(io)
    println(io, "source_dir: ", abspath(source_dir))
    println(io, "policy_file: ", abspath(policy_file))
    println(io, "functions_found: ", length(infos))
    println(io, "edges_found: ", sum(length(info.calls) for info in infos))
    println(io)
    println(io, "Modules discovered:")
    for m in sort(unique(module_of(info.full_name) for info in infos))
        println(io, "  - ", m)
    end
    println(io)

    if isempty(violations)
        println(io, "STATUS: PASS")
        println(io)
        println(io, "No architecture policy violations were found.")
    else
        println(io, "STATUS: FAIL")
        println(io)
        println(io, "Violations found: ", length(violations))
        println(io)
        for (i, v) in enumerate(violations)
            println(io, "[$i] ", v.rule)
            println(io, "    source: ", v.src)
            println(io, "    target: ", v.dst)
            println(io, "    detail: ", v.detail)
            println(io)
        end
    end
end

println(read(report_file, String))

isempty(violations) || exit(1)
