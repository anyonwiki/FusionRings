#!/usr/bin/env julia

using FuncDeps
using TOML

function usage()
  return println(
    """
Usage:
  julia --project=ci/FuncDepsEnv ci/run_funcdeps.jl <source_dir> <output_dir> [policy_file]

Example:
  julia --project=ci/FuncDepsEnv ci/run_funcdeps.jl src artifacts/funcdeps ci/architecture_policy.toml
""",
  )
end

function render_svg(cmd::Cmd, description::String)
  try
    run(cmd)
    println("Rendered ", description)
  catch err
    @warn "Could not render $description" exception=(err, catch_backtrace())
  end
end

function maybe_policy(path::Union{Nothing, String})
  path === nothing && return Dict{String, Any}()
  isfile(path) || return Dict{String, Any}()
  return TOML.parsefile(path)
end

function string_list(x)
  x === nothing && return String[]
  x isa AbstractVector || return String[]
  return String.(x)
end

length(ARGS) >= 2 || (usage(); exit(2))

source_dir = ARGS[1]
output_dir = ARGS[2]
policy_file = length(ARGS) >= 3 ? ARGS[3] : nothing

isdir(source_dir) || error("source_dir does not exist: $source_dir")
mkpath(output_dir)

policy = maybe_policy(policy_file)

println("Scanning source directory: ", abspath(source_dir))
infos = scan_project(source_dir)
idx = build_index(infos)

full_dot   = joinpath(output_dir, "funcdeps_full.dot")
cross_dot  = joinpath(output_dir, "funcdeps_cross_module.dot")
module_dot = joinpath(output_dir, "moduledeps.dot")
html_out   = joinpath(output_dir, "interactive_graph.html")

write_full_dot(infos, full_dot; cluster_by_module = false)
write_cross_module_dot(infos, cross_dot; cluster_by_module = true)
write_module_dot(infos, module_dot)
write_interactive_html(infos, html_out)

println("Wrote DOT/HTML outputs:")
println("  ", abspath(full_dot))
println("  ", abspath(cross_dot))
println("  ", abspath(module_dot))
println("  ", abspath(html_out))

# Focus views from policy.
focus_dir = joinpath(output_dir, "focus")
mkpath(focus_dir)

for mod in string_list(get(policy, "independent_modules", String[]))
  out = joinpath(focus_dir, "module_focus_$(replace(mod, r"[^A-Za-z0-9_]+" => "_")).dot")
  write_module_focus_dot(infos, mod, out)
  println("  ", abspath(out))
end

for file in string_list(get(policy, "independent_files", String[]))
  out = joinpath(focus_dir, "file_focus_$(replace(file, r"[^A-Za-z0-9_]+" => "_")).dot")
  write_file_focus_dot(infos, file, out)
  println("  ", abspath(out))
end

ind_funcs = get(policy, "independent_functions", Dict{String, Any}())
if ind_funcs isa Dict
  for f in vcat(
    string_list(get(ind_funcs, "cross_module", String[])),
    string_list(get(ind_funcs, "no_internal_calls", String[])),
  )
    out = joinpath(focus_dir, "function_focus_$(replace(f, r"[^A-Za-z0-9_]+" => "_")).dot")
    write_function_focus_dot(infos, f, out; depth = 2)
    println("  ", abspath(out))
  end
end

# Render SVGs when Graphviz is available.
dot = Sys.which("dot")
sfdp = Sys.which("sfdp")

if dot !== nothing
  render_svg(
    `$dot -Tsvg $module_dot -o $(joinpath(output_dir, "moduledeps.svg"))`, "moduledeps.svg"
  )
  render_svg(
    `$dot -Tsvg $cross_dot -o $(joinpath(output_dir, "funcdeps_cross_module.svg"))`,
    "funcdeps_cross_module.svg",
  )
else
  @warn "Graphviz dot not found; skipping module/cross-module SVG rendering"
end

if sfdp !== nothing
  render_svg(
    `$sfdp -Goverlap=prism -Gsplines=line -Gpack=true -Grepulsiveforce=2 -Tsvg $full_dot -o $(joinpath(output_dir, "funcdeps_full.svg"))`,
    "funcdeps_full.svg",
  )
elseif dot !== nothing
  render_svg(
    `$dot -Tsvg $full_dot -o $(joinpath(output_dir, "funcdeps_full.svg"))`,
    "funcdeps_full.svg via dot fallback",
  )
else
  @warn "Graphviz sfdp/dot not found; skipping full graph SVG rendering"
end

if dot !== nothing && isdir(focus_dir)
  for f in readdir(focus_dir; join = true)
    endswith(f, ".dot") || continue
    out = replace(f, r"\.dot$" => ".svg")
    render_svg(`$dot -Tsvg $f -o $out`, basename(out))
  end
end

summary = joinpath(output_dir, "funcdeps_summary.txt")
open(summary, "w") do io
  println(io, "FuncDeps summary")
  println(io, "================")
  println(io)
  println(io, "source_dir: ", abspath(source_dir))
  println(io, "functions_found: ", length(infos))
  println(io, "edges_found: ", sum(length(info.calls) for info in infos))
  println(io)
  println(io, "Modules:")
  for m in sort(unique(module_of(info.full_name) for info in infos))
    println(io, "  - ", m)
  end
  println(io)
  println(io, "Sample cross-module edges:")
  shown = 0
  for info in infos
    for callee in cross_module_callees_of(idx, info.full_name)
      println(io, "  ", info.full_name, " -> ", callee)
      shown += 1
      shown >= 30 && break
    end
    shown >= 30 && break
  end
  return shown == 0 && println(io, "  (none)")
end

println("Summary written to: ", abspath(summary))
