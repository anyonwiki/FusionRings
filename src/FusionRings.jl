module FusionRings

using Oscar, JSON, Base.Threads, Accessors
using LinearAlgebra: eigen, eigvals, diag
using ProgressMeter
using UUIDs
using Pkg
using Pkg.Artifacts

import Oscar:
  multiplication_table,
  is_commutative,
  rank,
  multiplicity,
  group,
  upper_central_series,
  is_nilpotent,
  tensor_product,
  automorphism_group,
  is_simple,
  is_integral
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
  !haskey(frd, [r, m, nnsd]) && error("code not found in database")
  !haskey(frd[[r, m, nnsd]], i) && error("code not found in database")

  return frd[[r, m, nnsd]][i]
end

export from_uuid

function from_uuid(s::String)
  !haskey(uuid_dict, s) && error("uuid not found")

  return uuid_dict[s]
end

function from_anyonwiki_code(v::AbstractVector{<:Integer})
  _ensure_frd_initialized()
  length(v) == 4 || error("anyonwiki_code expects a vector of 4 integers.")
  r, m, nnsd, i = Int.(collect(v))
  return from_anyonwiki_code(r, m, nnsd, i)
end

const fawc = from_anyonwiki_code
const awc  = anyonwiki_code

function __init__()
  # GLOBAL VARIABLES
  global QQb     = algebraic_closure(QQ)
  global QQab, ζ = abelian_closure(QQ)


  fusionrings_data_path       = joinpath( artifact"FusionRings", "FusionRings" )
  algebraic_numbers_data_path = joinpath( artifact"AlgebraicNumbers", "AlgebraicNumbers" )

  # IMPORT FUSION RINGS
  is_fr_data_file( s::String ) = s[1:4] == "mult" && s[end-4:end] == ".json"
  fr_filenames = filter( is_fr_data_file, readdir(fusionrings_data_path) )

  global fusion_ring_list =
    vcat([ import_rings(joinpath(fusionrings_data_path, fn)) for fn in fr_filenames ]...);

  global frl = fusion_ring_list

  global uuid_dict = Dict(uuid(r) => r for r in frl)

  # for unknown rings, the rank, mult, and nnsd can
  # be determined quickly. We will group the known fusion rings by
  # those properties so its faster to identify unknown rings.
  grouped_by_first3 = Dict{Vector{Int64}, Vector{FusionRing}}()
  for ring in frl
    key = (ring.anyonwiki_code)[1:3]
    if !haskey(grouped_by_first3, key)
      grouped_by_first3[key] = FusionRing[]
    end
    push!(grouped_by_first3[key], ring)
  end

  fourth_to_dict(v::Vector{FusionRing}) = Dict((r.anyonwiki_code)[4] => r for r in v)

  global fusion_ring_dict = Dict(k => fourth_to_dict(v) for (k, v) in grouped_by_first3)

  global frd = fusion_ring_dict

  # CHECK IF FIRST TIME USING PACKAGE

  #split number data in separate files
  # these will be loaded on demand rather than all at the same time.
  global datadir = joinpath(@__DIR__, "data")
  fns = readdir(datadir)

  ids = Oscar.load(joinpath(algebraic_numbers_data_path, "qqb_ids.mrdi"))
  nf = length(ids)
  if "split_number_data" ∉ fns || length(readdir(joinpath(datadir, "split_number_data"))) < nf
    println("Dataset of algebraic numbers not yet optimized. Optimizing for future use.")

    # Create directory for numbers
    splitdatapath = joinpath(datadir, "split_number_data/")
    mkpath(splitdatapath)

    println("Importing algebraic numbers.")
    # Import qqb numbers from big files
    qqb_nums = begin
      nums = Oscar.load(joinpath(algebraic_numbers_data_path, "qqb_vals.mrdi"))

      [(ids[i], nums[i]) for i in 1:length(ids)]
    end

    println()
    # Export qqb numbers
    function exportnum(tuple)
      path = joinpath(datadir, "split_number_data/", tuple[1]*".mrdi")
      return Oscar.save(path, tuple[2])
    end

    @showprogress desc="Exporting numbers separately." dt = 0.5 for num in qqb_nums
      exportnum(num)
    end
    println("Dataset is optimized.")
  end

  # qqb_dict starts empty and grows on demand
  global qqb_dict = Dict{String, QQBarFieldElem}()

end

end
