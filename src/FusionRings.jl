module FusionRings

using Oscar, JSON, Base.Threads, Accessors
using LinearAlgebra: eigen, eigvals, diag
import Oscar:
  multiplication_table, is_commutative, rank, multiplicity, group, upper_central_series
import Base: sort, names

include("general_functions.jl")
include("structs.jl")
include("formatting_and_printing.jl")
include("creation.jl")
include("properties.jl")
include("operations.jl")
include("import_data.jl")

export qqb_dict, fusion_ring_list, frl, fusion_ring_dict, frd, from_anyonwiki_code, fawc

"""Internal helper: ensure ring lookup dictionary has been initialized."""
function _ensure_frd_initialized()
  isdefined(@__MODULE__, :frd) ||
    error("FusionRings data not initialized (frd is undefined).")
  frd isa Dict || error("frd is defined but not a Dict).")
  return nothing
end

function is_optimized_db()
  datadir = joinpath(@__DIR__, "data")
  return ""
end

"""
    from_anyonwiki_code(r, m, nnsd, i) -> FusionRing
    from_anyonwiki_code(v::AbstractVector{<:Integer}) -> FusionRing

Lookup a stored fusion ring by its AnyonWiki code `[r, m, nnsd, i]`.
This requires package data to be loaded (normally happens during `__init__`).
"""
function from_anyonwiki_code(r::Integer, m::Integer, nnsd::Integer, i::Integer)
  _ensure_frd_initialized()
  return frd[Int[r, m, nnsd]][i]
end

function from_anyonwiki_code(v::AbstractVector{<:Integer})
  _ensure_frd_initialized()
  length(v) == 4 || error("anyonwiki_code expects a vector of 4 integers.")
  r, m, nnsd, i = Int.collect(v)
  return frd[r, m, nnsd][i]
end

const fawc = from_anyonwiki_code
const awc  = anyonwiki_code

function __init__()
  # GLOBAL VARIABLES
  global QQb     = algebraic_closure(QQ)
  global QQab, ζ = abelian_closure(QQ)

  global datadir = joinpath(@__DIR__, "data")
  fns = readdir(datadir)

  # IF FIRST TIME USING PACKAGE: split number data in separate files
  # these will be loaded on demand rather than all at the same time.

  if "split_number_data" ∉ fns
    println("Dataset of algebraic numbers not yet optimized. Optimizing for future use.")
    # Create directory for numbers
    splitdatapath = joinpath(datadir, "split_number_data/")
    mkdir(splitdatapath)

    println("Importing algebraic numbers.")
    # Import qqb numbers from big files
    qqb_nums = begin
      ids  = Oscar.load(joinpath(datadir, "qqb_ids.mrdi"))
      nums = Oscar.load(joinpath(datadir, "qqb_vals.mrdi"))

      [ids[i], nums[i] for i in 1:length(ids)]
    end

    println("Exporting numbers separately.")
    # Export qqb numbers
    function exportnum(tuple)
      dir = joinpath(datadir, "split_number_data/", tuple[1]*".mrdi")
      return Oscar.save(path, tuple[2])
    end

    exportnum.(qqb_nums)
    println("Dataset is optimized.")
  end

  # qqb_dict starts empty and grows on demand
  global qqb_dict = Dict{String, QQBarFieldElem}()

  # IMPORT FUSION RINGS
  # TODO: how long does sorting take? Might want to store the permutation vector 
  # once and reuse it
  global fusion_ring_list = sort( # Stored list is unsorted so we still need to sort
    import_rings(joinpath(datadir, "fusionrings.json"));
    by = (x -> (x.anyonwiki_code)[[2, 1, 3, 4]]),
  )
  global frl = fusion_ring_list

  # for unknown rings the first 3 indices of the anyonwiki_code can 
  # be determined quickly. We will group the known fusion rings by 
  # the first 3 indices and then, separately, by the 4th

  grouped_by_first3 = Dict{Vector{Int64}, Vector{FusionRing}}()
  for ring in frl
    key = anyonwiki_code(ring)[1:3]
    if !haskey(grouped_by_first3, key)
      grouped_by_first3[key] = FusionRing[]
    end
    push!(grouped_by_first3[key], ring)
  end

  fourth_to_dict(v::Vector{FusionRing}) = Dict((r.anyonwiki_code)[4] => r for r in v)

  global fusion_ring_dict = Dict(k => fourth_to_dict(v) for (k, v) in grouped_by_first3)
  #Dict( anyonwiki_code(r) => r for r in frl )
  return global frd = fusion_ring_dict
end

end
