#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                   EXPORTING AND IMPORTING QQBARFIELDELEMS                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# The QQBarFieldElem objects are quite heavy to load so we will load them all in a dictionary 
# and provide functions to convert qqbar elems to keys and vice versa

# Generate unique ID for a QQBarFieldElem
export qqb_id

function qqb_id(x::QQBarFieldElem)
  mp = minimal_polynomial(x)
  coeffs = string.(collect(coefficients(mp)))
  us = fill("_", degree(mp) + 1)

  numstring = string(rootnum(x))

  return stringriffle(coeffs, us) * "_" * numstring
end

function qqb_id(arr::Array{QQBarFieldElem})
  return qqb_id.(arr)
end

function qqb_id(melem::AbstractAlgebra.Generic.MatSpaceElem{QQBarFieldElem})
  return qqb_id.(Matrix(melem))
end

function rootnum(x::QQBarFieldElem)::Int64
  p   = minimal_polynomial(x)
  rts = roots(QQBar, p)
  sr  = sort(rts; by = root_sort_crit)
  return findfirst(y -> y == x, sr)
end

# This is the sort criterion for roots used by mathematica 
# and by the anyonwiki on 28/12/2025
function root_sort_crit(x)
  return (- Int(is_real(x)), real(x), imag(x))
end

function save_qqb_num(
  x::QQBarFieldElem;
  directory::AbstractString = joinpath(datadir, "split_number_data"),
)
  id = qqb_id(x)
  mkpath(directory)
  Oscar.save(joinpath(directory, id*".mrdi"), x)
  qqb_dict[id] = x
  return id
end

# Loading from library of qqb elements
# If number already loaded, use that one, otherwise load and
# add to qqb_dict
# If id not available, create qqb num and save it

export from_qqb_id

function from_qqb_id(s::String)
  if haskey(qqb_dict, s)
    return qqb_dict[s]
  else
    fn = joinpath(datadir, "split_number_data", s*".mrdi")

    if isfile(fn)
      val = Oscar.load(fn)
      qqb_dict[s] = val
      return val
    else
      spl = split(s, "_")
      l   = length(spl)
      n   = parse(Int64, last(spl))
      cfs = ZZ.(parse.(Int64, spl[1:(end - 2)]))

      qqbval = sort(roots(QQBar, polynomial(ZZ, cfs)); by = root_sort_crit)[n]

      save_qqb_num(qqbval)

      return qqbval
    end
  end
end

function from_qqb_id(a::Vector)
  return from_qqb_id.(a)
end

function from_qqb_id(a::Array)
  return from_qqb_id.(a)
end

function from_qqb_id(a::Matrix)
  return from_qqb_id.(a)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           IMPORTING FUSION RINGS                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
# The fusion rings are stored as json files. Not all data
# types (e.g. complex numbers) are supported by JSON so we 
# store those using a variety of hacks.
# The following functions convert the stored data back 
# to their proper types.
#
# TODO: some of the if clauses below are necessary for legacy
# compatibility. Once all json files have the correct format
# we should remove it since it slows down the import

# helper function to deal with missing keys and values equal to nothing
function safe_fetch(js, key)
  ks = keys(js)

  key ∉ ks && return missing

  val = js[key]

  isnothing(val) && return missing

  return val
end

# mult tab
function mtfromjs(js)::Array{Int64, 3}
  jsmt = safe_fetch(js, "mult_tab")
  ismissing(jsmt) && throw(ArgumentError("mult_tab is required"))
  jsmt isa AbstractVector || throw(ArgumentError("mult_tab must be a cubic array"))

  r = length(jsmt)
  r > 0 || throw(ArgumentError("mult_tab must have positive rank"))
  all(
    plane ->
      plane isa AbstractVector &&
      length(plane) == r &&
      all(row -> row isa AbstractVector && length(row) == r, plane),
    jsmt,
  ) || throw(ArgumentError("mult_tab must be a cubic array"))

  mt = zeros(Int, r, r, r)
  for i in 1:r, j in 1:r, k in 1:r
    value = jsmt[i][j][k]
    value isa Integer ||
      throw(ArgumentError("every mult_tab entry must be an integer"))
    mt[i, j, k] = Int(value)
  end
  return mt
end

function uuidfromjs(js)::Union{String, Missing}
  return safe_fetch(js, "uuid")
end

# anyonwiki_code
function fcfromjs(js)::Union{Missing, Vector{Int64}}
  k = keys(js)

  if "anyonwiki_id" ∈ k
    id = js["anyonwiki_id"]

    isnothing(id) && return missing

    return id_to_formal_code(id)
  end

  fc = if "formal_code" ∈ k
    js["formal_code"]
  elseif "anyonwiki_code" ∈ k
    js["anyonwiki_code"]
  else
    nothing
  end

  if isnothing(fc) || length(fc) == 0
    return missing
  else
    return fc
  end
end

function formal_code_to_id(code::Vector{Int64})
  return "anyonwiki_fcrm_fr__" * join(string.(code), "_")
end

function id_to_formal_code(id::String)
  c = parse.(Int64, split(id[20:end], "_"))

  length(c) ≠ 4 && error("Invalid fusion ring id")

  return c
end

# importing names
_nonames = Dict(
  "quantum_group_like" => missing,
  "group_like"         => missing,
  "physics"            => missing,
  "miscellaneous"      => missing,
)

function mscnames(v::Vector)
  return Dict(
    "quantum_group_like" => missing,
    "group_like"         => missing,
    "physics"            => missing,
    "miscellaneous"      => string.(v),
  )
end

# names
function nfromjs(js)
  nms = safe_fetch(js, "names")

  ismissing(nms) && return _nonames

  nms isa Vector && isempty(nms) && return _nonames

  nms isa Vector && return mscnames(nms)

  return nms
end

# texnames
function tnfromjs(js)
  nms = safe_fetch(js, "texnames")

  ismissing(nms) && return _nonames

  nms isa Vector && isempty(nms) && return _nonames

  nms isa Vector && return mscnames(nms)

  return nms
end

# labels
function lfromjs(js)
  lbs = safe_fetch(js, "labels")

  ismissing(lbs) && return missing
  lbs isa AbstractVector || throw(ArgumentError("labels must be an array or null"))

  return string.(lbs)
end

# characters
function chfromjs(js)
  chrs = safe_fetch(js, "characters")
  ismissing(chrs) && return missing
  chrs isa AbstractVector ||
    throw(ArgumentError("characters must be a square array or null"))

  r = length(chrs)
  r > 0 || throw(ArgumentError("characters cannot be an empty array"))
  all(row -> row isa AbstractVector && length(row) == r, chrs) ||
    throw(ArgumentError("characters must be a square array"))

  return [string(chrs[i][j]) for i in 1:r, j in 1:r]
end

# sub-fusion rings
function sfrfromjs(js)
  sfr = safe_fetch(js, "non_trivial_sub_fusion_rings")
  ismissing(sfr) && return missing
  sfr isa AbstractVector ||
    throw(ArgumentError("non_trivial_sub_fusion_rings must be an array or null"))

  decoded = Tuple{Vector{Int64}, String}[]
  for entry in sfr
    entry isa AbstractVector && length(entry) == 2 ||
      throw(ArgumentError("each stored subring must be [indices, uuid]"))
    indices, id = entry
    indices isa AbstractVector && all(x -> x isa Integer, indices) ||
      throw(ArgumentError("stored subring indices must be integers"))
    id isa AbstractString ||
      throw(ArgumentError("stored subring UUID must be a string"))
    push!(decoded, (Int64.(indices), String(id)))
  end
  return decoded
end

function fpdfromjs(js)
  fpd = safe_fetch(js, "frobenius_perron_dimension")
  ismissing(fpd) && return missing
  fpd isa AbstractString ||
    throw(ArgumentError("frobenius_perron_dimension must be a QQBar ID string or null"))
  return String(fpd)
end

function fpdsfromjs(js)
  fpds = safe_fetch(js, "frobenius_perron_dimensions")
  ismissing(fpds) && return missing
  fpds isa AbstractVector || throw(
    ArgumentError("frobenius_perron_dimensions must be an array of QQBar ID strings or null")
  )
  all(x -> x isa AbstractString, fpds) ||
    throw(ArgumentError("every stored FP dimension must be a QQBar ID string"))

  return String.(fpds)
end

function fcdfromjs(js)
  fcd = safe_fetch(js, "formal_codegrees")
  ismissing(fcd) && return missing
  fcd isa AbstractVector ||
    throw(ArgumentError("formal_codegrees must be an array of QQBar ID strings or null"))
  all(x -> x isa AbstractString, fcd) ||
    throw(ArgumentError("every formal codegree must be a QQBar ID string"))

  return String.(fcd)
end

function ctpfromjs(js)
  hcwp = safe_fetch(js, "has_categories_with_props")

  if hcwp isa Vector
    return Dict{String, Vector{Any}}(
      row[1] => [row[2], row[3][1], row[3][2]] for row in hcwp
    )
  end

  if ismissing(hcwp)
    return Dict(
      "Fusion"    => [missing, "", "Unknown to AnyonWiki"],
      "Pivotal"   => [missing, "", "Unknown to AnyonWiki"],
      "Unitary"   => [missing, "", "Unknown to AnyonWiki"],
      "Spherical" => [missing, "", "Unknown to AnyonWiki"],
      "Braided"   => [missing, "", "Unknown to AnyonWiki"],
      "Ribbon"    => [missing, "", "Unknown to AnyonWiki"],
      "Modular"   => [missing, "", "Unknown to AnyonWiki"],
    )
  end

  return hcwp
end

function ctsfromjs(js)
  return safe_fetch(js,"categorifications")
end

function refsfromjs(js)
  refs = safe_fetch(js, "references")

  refs isa Vector && return Dict{String, Vector{String}}("All" => refs)

  return refs
end

function sfromjs(js)
  sftw = safe_fetch(js, "software")

  sftw isa Vector && return Dict{String, Vector{String}}("All" => sftw)

  return sftw
end

function agfromjs(js)
  ag = safe_fetch(js, "all_gradings")
  ismissing(ag) && return missing
  ag isa AbstractVector ||
    throw(ArgumentError("all_gradings must be an array or null"))

  decoded = Tuple{Vector{Int64}, String}[]
  for grading in ag
    grading isa AbstractVector && length(grading) == 2 ||
      throw(ArgumentError("each stored grading must be [degrees, uuid]"))
    degrees, id = grading
    degrees isa AbstractVector && all(x -> x isa Integer, degrees) ||
      throw(ArgumentError("stored grading degrees must be integers"))
    id isa AbstractString ||
      throw(ArgumentError("stored grading-ring UUID must be a string"))
    push!(decoded, (Int64.(degrees), String(id)))
  end
  return decoded
end

function ucsfromjs(js)
  ucs = safe_fetch(js, "upper_central_series")
  ismissing(ucs) && return missing
  ucs isa AbstractVector ||
    throw(ArgumentError("upper_central_series must be an array or null"))

  decoded = Tuple{Vector{Int64}, String}[]
  for term in ucs
    term isa AbstractVector && length(term) == 2 ||
      throw(ArgumentError("each upper-central-series term must be [indices, uuid]"))
    indices, id = term
    indices isa AbstractVector && all(x -> x isa Integer, indices) ||
      throw(ArgumentError("upper-central-series indices must be integers"))
    id isa AbstractString ||
      throw(ArgumentError("upper-central-series UUID must be a string"))
    push!(decoded, (Int64.(indices), String(id)))
  end
  return decoded
end

function rfromjs(js)
  stored = safe_fetch(js, "realizations")
  ismissing(stored) && return Dict{String, Any}("tensor_product" => missing)
  stored isa AbstractDict ||
    throw(ArgumentError("realizations must be a JSON object or null"))

  result = Dict{String, Any}(String(k) => v for (k, v) in stored)
  if !haskey(result, "tensor_product") && haskey(result, "TensorProduct")
    result["tensor_product"] = pop!(result, "TensorProduct")
  end
  get!(result, "tensor_product", missing)
  return result
end

function afromjs(js)
  aut = safe_fetch(js, "automorphisms")
  ismissing(aut) && return missing
  aut isa AbstractDict ||
    throw(ArgumentError("automorphisms must be a JSON object or null"))
  haskey(aut, "cycles") ||
    throw(ArgumentError("automorphisms must contain a cycles field"))

  generators = aut["cycles"]
  generators isa AbstractVector ||
    throw(ArgumentError("automorphism cycles must be an array"))

  mt = safe_fetch(js, "mult_tab")
  ismissing(mt) && throw(ArgumentError("mult_tab is required to decode automorphisms"))
  rk = length(mt)
  isempty(generators) &&
    return permutation_group(rk, [perm(collect(1:rk))])

  converted = map(generators) do generator
    generator isa AbstractVector ||
      throw(ArgumentError("each automorphism generator must be an array of cycles"))
    all(
      cycle -> cycle isa AbstractVector && all(x -> x isa Integer, cycle), generator
    ) || throw(ArgumentError("automorphism cycles must contain integer indices"))
    return cperm([Int64.(cycle) for cycle in generator])
  end

  return permutation_group(rk, converted)
end

export import_ring

function import_ring(filename::String; skip_check = false)
  return import_ring(JSON.parsefile(filename); skip_check = skip_check)
end

function import_ring(js::JSON.Object{String, Any}; skip_check = false)
  mt = mtfromjs(js)
  rk = size(mt, 1)

  return fusion_ring(
    mt;
    uuid                        = uuidfromjs(js),
    anyonwiki_code              = fcfromjs(js),
    names                       = nfromjs(js),
    texnames                    = tnfromjs(js),
    labels                      = lfromjs(js),
    characters                  = chfromjs(js),
    sub_fusion_rings            = sfrfromjs(js),
    frobenius_perron_dimension  = fpdfromjs(js),
    frobenius_perron_dimensions = fpdsfromjs(js),
    formal_codegrees            = fcdfromjs(js),
    has_categories_with_props   = ctpfromjs(js),
    categorifications           = ctsfromjs(js),
    references                  = refsfromjs(js),
    software                    = sfromjs(js),
    all_gradings                = agfromjs(js),
    upper_central_series        = ucsfromjs(js),
    realizations                = rfromjs(js),
    automorphism_group          = afromjs(js),
    skip_check                  = skip_check,
  )
end

export import_rings

function import_rings(filename::String; skip_check = false)
  jsdict = JSON.parsefile(filename,Dict{String,Any});

  frlist = FusionRing[]

  k = keys(jsdict)
  indices = "order" ∈ k ? jsdict["order"] : eachindex(jsdict["data"])

  @showprogress dt=1 desc = "Importing fusion rings" for ind in indices
    js = jsdict["data"][ind]
    push!(frlist, import_ring(js, skip_check = skip_check))
  end

  return frlist
end

############################################################
# Exporting Fusion rings
############################################################

function missing_to_nothing(x)
  !ismissing(x) && return x

  return nothing
end

function mttojs(fr::FusionRing)::Vector{Vector{Vector{Int64}}}
  mt = multiplication_table(fr)
  r  = rank(fr)
  return [[[mt[i, j, k] for k in 1:r] for j in 1:r] for i in 1:r]
end

function actojs(fr::FusionRing)::Union{Vector{Int64}, Nothing}
  c = anyonwiki_code(fr)
  return ismissing(c) ? nothing : c
end

function chtojs(fr::FusionRing)
  ch = fr.characters
  if ch !== missing
    r = rank(fr)
    return [[ch[i, j] for j in 1:r] for i in 1:r]
  else
    return nothing
  end
end

function sfrtojs(fr::FusionRing)
  sfr = fr.sub_fusion_rings
  (ismissing(sfr) || isnothing(sfr)) && return nothing

  return [[t[1], t[2]] for t in sfr]
end

function fpdtojs(fr::FusionRing)::Union{String, Nothing}
  fpd = fr.frobenius_perron_dimension
  (ismissing(fpd) || isnothing(fpd)) && return nothing

  return fpd
end

function fpdstojs(fr::FusionRing)::Union{Vector{String}, Nothing}
  fpd = fr.frobenius_perron_dimensions
  (ismissing(fpd) || isnothing(fpd)) && return nothing

  return fpd
end

function fcdtojs(fr::FusionRing)::Union{Vector{String}, Nothing}
  fcd = fr.formal_codegrees
  (ismissing(fcd) || isnothing(fcd)) && return nothing

  return fcd
end

function nfcdtojs(fr::FusionRing)::Union{Vector{Vector{Float64}}, Nothing}
  fcd = fr.formal_codegrees
  (ismissing(fcd) || isnothing(fcd)) && return nothing

  return reim.(ComplexF64.(from_qqb_id(fcd)))
end

function nchtojs(fr::FusionRing)::Union{Vector{Vector{Vector{Float64}}}, Nothing}
  nc = fr.characters
  (ismissing(nc) || isnothing(nc)) && return nothing

  r = rank(fr)
  splitchars = reim.(numeric_characters(fr))
  return [[splitchars[i, j] for j in 1:r] for i in 1:r]
end

function nfpdtojs(fr::FusionRing)
  fpd = fr.frobenius_perron_dimension
  (ismissing(fpd) || isnothing(fpd)) && return nothing

  return reim(numeric_fpdim(fr))
end

function nfpdstojs(fr::FusionRing)
  fpd = fr.frobenius_perron_dimensions
  (ismissing(fpd) || isnothing(fpd)) && return nothing

  return reim.(numeric_fpdims(fr))
end

function reim(x::ComplexF64)::Vector{Float64}
  return [real(x), imag(x)]
end

function reim(x::Float64)::Vector{Float64}
  return [x, 0.0]
end

function reim(mat::Matrix{ComplexF64})::Vector{Vector{Vector{Float64}}}
  return [[reim(r[i]) for i in eachindex(r)] for r in eachrow(mat)]
end

function reim(vv::Vector{Vector{ComplexF64}})::Vector{Vector{Vector{Float64}}}
  return [[reim(coef) for coef in vec] for vec in vv]
end

#TODO: for some rings we automatically know these props e.g.
# abelian groups: all true (exept maybe modular?)
# nonabelian groups: fusion, piv, unitary, spherical true, rest false
# quantum group like rings: need to look this up but a lot is known as well
#

function ctstojs(fr::FusionRing)
  return missing_to_nothing(fr.categorifications)
end

function agtojs(fr::FusionRing)
  ag = fr.all_gradings
  (ismissing(ag) || isnothing(ag)) && return nothing

  gap_id(id::String) = small_group_identification(to_group(from_uuid(id)))

  #TODO: add gap_id to info

  return [[g[1], g[2]] for g in ag]
end

function ucstojs(fr::FusionRing)
  ucs = fr.upper_central_series
  (ismissing(ucs) || isnothing(ucs)) && return nothing

  return ucs
end

function rtojs(fr::FusionRing)
  function convert_realization(val)
    ismissing(val) && return nothing
    return [[uuid(ring) for ring in list] for list in val]
  end

  return Dict(k => convert_realization(v) for (k, v) in realizations(fr))
end

function auttojs(fr::FusionRing)
  ag = automorphism_group(fr)

  (ismissing(ag) || isnothing(ag)) && return nothing

  mgs = minimal_size_generating_set(ag)
  remove_unit_cycles(cyc) = filter(x -> length(x) != 1, cyc)
  cyc = remove_unit_cycles.(cycles.(mgs))

  gap_id = small_group_identification(ag)

  return Dict("group" => tex_describe(ag), "cycles" => cyc)
end

function tex_describe(gp::Group)
  s = describe(gp)
  s == "1" && return "Trivial"
  # Order matters: replace longer/specific patterns first
  s = replace(s, "QDn" => "QD_{n}")
  s = replace(s, "Q8" => "\\mathbb{H}")
  s = replace(s, r"S\d+" => m -> "S_{$(m[2:end])}")
  s = replace(s, r"Q\d+" => m -> "Q_{$(m[2:end])}")
  s = replace(s, "Z" => "\\mathbb{Z}")
  s = replace(s, r"C\d+" => m -> "\\mathbb{Z}_{$(m[2:end])}")
  s = replace(s, r"D\d+" => m -> "D_{$(m[2:end])}")
  s = replace(s, "x" => "\\times")
  return s
end

function intojs(fr::FusionRing)::Union{Bool, Nothing}
  ucs = fr.upper_central_series
  (ismissing(ucs) || isnothing(ucs)) && return nothing

  return is_nilpotent(fr)
end

function istojs(fr::FusionRing)::Union{Bool, Nothing}
  sfr = fr.sub_fusion_rings
  (ismissing(sfr) || isnothing(sfr)) && return nothing

  return is_simple(fr)
end

function iitojs(fr::FusionRing)::Union{Bool, Nothing}
  fpds = fr.frobenius_perron_dimensions

  (ismissing(fpds) || isnothing(fpds)) && return nothing

  return is_integral(fr)
end

function iwitojs(fr::FusionRing)::Union{Bool, Nothing}
  fpd = fr.frobenius_perron_dimension

  (ismissing(fpd) || isnothing(fpd)) && return nothing

  return is_weakly_integral(fr)
end

function intgtojs(fr::FusionRing)::Union{Bool, Nothing}
  ag = fr.all_gradings

  (ismissing(ag) || isnothing(ag)) && return nothing

  return length(ag) > 1
end

function write_json(filename::String, data::Dict)
  open(filename, "w") do f
    return JSON.json(f, data; pretty = true, inline_limit = 10)
  end
end

ringfieldinfo = """
    The following conventions are used in explaining the values of the dictionary:
    * A qqb_id is a string that uniquely describes an algebraic number. It is formatted as a list of n integers a0, ..., an, separated by underscores, followed by a double underscore, followed by an integer i: "a0_a1_..._an__i". Here a0 to an are the coefficients of a polynomial a0 + a1*x + ... + an*x^n and i denotes the i'th root of that polynomial. The indexing of roots of the polynomial takes the real roots first, in increasing order. Then come complex conjugate pairs of roots, sorted first by increasing real part and second by increasing complex part.
    * A ctf, or complex tuple of floats, is a representation of a complex floating point number a + i b by vector [ a, b ] where a and b are real floating point numbers.
    * If A is a JSON array then A[i] is its i'th element.

    * The parent ring is the current ring being represented by the JSON dictionary. When talking about subrings and gradings, the parent ring represents the ring of which we're interested in its subrings and gradings.
    * Every fusion ring is uniquely identified by a uuid. By ring[ring_uuid] we mean the fusion ring with uuid equal to ring_uuid
    * All stored rings have a fixed order of their elements. Therefore every element of a fusion ring can be represented by a positive integer from 1 to r = rank(ring) we will often use elements and indices representing those elements interchangeably. By the i'th element of the ring[uuid] we mean the i'th element for the stored ring unless stated otherwise.
    * A value of null for any field means the data is missing. In particular, it does not imply that the data doesn't exist.
    * For all technical definitions we refer to doi: 10.1090/surv/205.

    The interpretation of the values of the fields of a fusion ring is the following.
    * mult_tab: triply nested array of integers representing the structure constants of the fusion ring: N_{a,b}^c = mult_tab[a][b][c]
    * uuid: UUID1 string that uniquely represents the fusion ring. It is independent of any property of the fusion ring and will therefore not change if a property is found to be incorrect, e.g., due to an incorrect property.
    * anyonwiki_code: a list of four integers r, m, nnsd, i, that identify a unique fusion ring. Here
      * r is the rank of the ring
      * m is the multiplicity of the ring
      * nnsd is the number of non-self dual elements of the fusion ring
      * is an arbitrary integer that distinguished between rings with the same values of r, m, and nnsd.
    * names: a JSON dictionary mapping naming conventions to lists of strings of names given using that convention. The conventions at the moment are
      * "quantum_group_like": names associated to quantum groups at level k, such as "psu(2)_5", "so(7)_2"
      * "group_like": names associated to the theory of finite groups, such as "Z_2", "Rep(D_6)". Names associated to near-group or group theoretical fusion rings do not belong here but in miscelaneous.
      * "physics": names associated to applications in physics, such as "Fibonacci", "Ising", "Potts".
      * "miscellaneous": names not belonging to another of the above categories.
    * texnames: a JSON dictionary mapping naming conventions to lists of strings of names typeset in LaTeX given using that convention. The conventions are the same as for the names field.
    * labels: list of strings used to label the elements of the fusion ring. This is purely cosmetic and has no influence on any other properties. The default is a list of bold digits from 1 to the rank of the ring.
    * characters: vector of vectors [ v_1, ..., v_r ] where each v_i is a vector of r qqb_ids representing the image of the elements of the fusion ring under the i'th character.
    * non_trivial_sub_fusion_rings: list of vectors [ els, uuid ] where els are the elements of the parent ring that form a subring isomorphic to ring[uuid].
    * frobenius_perron_dimension: qqb_id of the Frobenius-Perron dimension of the fusion ring.
    * frobenius_perron_dimensions: vector of qqb_ids of the Frobenius-Perron dimensions of the elements of the fusion ring.
    * formal_codegrees: vector of qqb_ids that are the eigenvalues of the matrix \\sum_{i=1}^{rank(fr)}N_{i^*} N_{i} where N_{i} is the left regular representation of left-multiplication by the i'th basis element and i^* is the dual element of i.
    * numeric_formal_codegrees: vector of ctfs that are the eigenvalues of the matrix \\sum_{i=1}^{rank(fr)}N_{i^*} N_{i} where N_{i} is the left regular representation of left-multiplication by the i'th basis element and i^* is the dual element of i.
    * numeric_characters: vector of vectors [ v_1, ..., v_r ] where each v_i is a vector of r ctfs representing the image of the elements of the fusion ring under the i'th character.
    * numeric_frobenius_perron_dimension: ctf of the Frobenius-Perron dimension of the fusion ring.
    * numeric_frobenius_perron_dimensions: vector of ctfs of the Frobenius-Perron dimensions of the elements of the fusion ring.
    * has_categories_with_props: a JSON dictionary whose keys are
      * "Fusion"
      * "Pivotal"
      * "Unitary"
      * "Spherical"
      * "Braided"
      * "Ribbon"
      * "Modular"
      and whose fields are lists [ bool, method, reason ] where
      * bool: equals true if the ring is categorifiable to a category with the respective property, false if it is known it doesn't.
      * method: can be "Theory" if the information of bool is based on a theoretical result,  "Computer:symbolic" if it is based on a computer calculation using exact symbolic methods, or "Computer:numeric" if it is based on a computer calculation using inexact numeric methods.
      * reason: gives a more in-depth reason for why the value of bool is what it is. This could, e.g., be a reference to a theorem in a paper or a version of a software package used.
    * categorifications: a list of uuids of known fusion categories that categorify the fusion ring. It only contains uuids of categories of which the data is stored.
    * references: JSON dictionary mapping names of fields to a list of references to the paper that played a significant role in obtaining the data in the way it is represented here. Special field names are the same as for software. Only papers that have lead to the data as currently represented are included and thus no papers that represent theory that was not directly used, or ,e.g. , data in another format that was not used to obtain current data.
    * software: JSON dictionary mapping names of fields to a list of reference to software that played a significant role in obtain the data in the way it is represented here. Special field names are
      * "all": when all fields of the ring point to the same software
      * "all_other_data": when all other data, besides the data having specific references, points to the same software.
    * all_gradings: vector of vectors [ els, uuid ] where ring[uuid] the group ring that grades this fusion ring, and els are elements of ring[uuid] that grade the elements of the parent ring.
    * upper_central_series: list of vectors v_i =  [ els_i, uuid_i ] where ring[uuid_i] is the ring isomorphic to the adjoint fusion ring of the ring[uuid_{i-1}]. els_i are the elements of ring[uuid_{i-1}] that form its adoint fusion ring. v1 is by definition the couple of all elements of the parent ring and the parent ring itself. Each adjoint ring has its elements in the same order as the original ring and thus not necessarily in the order of the stored ring.
    * realizations: JSON dictionary mapping strings representing realizations of fusion rings in terms of other ones to data that allows to reconstruct the realization. At the moment it contains the following fields
      * "tensor_product": vector of vectors of uuids representing rings whose (based) tensor product is isomorphic to this ring.
    * automorphisms: JSON dictionary with the following fields:
      * "group": string containing a LaTeX name of the automorphism group of the fusion ring.
      * "cycles": minimal list of vectors of integers representing cycles that generate the automorhism group of the fusion ring with the current order of elements.
    * is_group_ring: true if the ring is a group ring, false if not.
    * is_nilpotent: true if the ring is nilpotent, false if not.
    * is_simple: true if the ring is simple, false if not.
    * is_integral: true if the ring is integral, false if not.
    * is_weakly_integral: true if the ring is nilpotent, false if not.
    * is_non_trivially_graded: true if the ring has a non-trivial grading, false if not.
    * is_commutative: true if the ring is commutative, false if not.
     """

function ring_to_dict(fr)
  infostring = "Fusion ring given by JSON dictionary.\n" * ringfieldinfo

  return Dict(
    "mult_tab"                            => mttojs(fr),
    "uuid"                                => uuid(fr),
    "anyonwiki_code"                      => actojs(fr),
    "names"                               => names(fr),
    "texnames"                            => tex_names(fr),
    "labels"                              => labels(fr),
    "characters"                          => chtojs(fr),
    "non_trivial_sub_fusion_rings"        => sfrtojs(fr),
    "frobenius_perron_dimension"          => fpdtojs(fr),
    "frobenius_perron_dimensions"         => fpdstojs(fr),
    "formal_codegrees"                    => fcdtojs(fr),
    "numeric_formal_codegrees"            => nfcdtojs(fr),
    "numeric_characters"                  => nchtojs(fr),
    "numeric_frobenius_perron_dimension"  => nfpdtojs(fr),
    "numeric_frobenius_perron_dimensions" => nfpdstojs(fr),
    "has_categories_with_props"           => categories_with_properties(fr),
    "categorifications"                   => ctstojs(fr),
    "references"                          => fr.references,
    "software"                            => fr.software,
    "all_gradings"                        => agtojs(fr),
    "upper_central_series"                => ucstojs(fr),
    "realizations"                        => rtojs(fr),
    "automorphisms"                       => auttojs(fr),
    "is_group_ring"                       => is_group_ring(fr),
    "is_nilpotent"                        => intojs(fr),
    "is_simple"                           => istojs(fr),
    "is_integral"                         => iitojs(fr),
    "is_weakly_integral"                  => iwitojs(fr),
    "is_non_trivially_graded"             => intgtojs(fr),
    "is_commutative"                      => is_commutative(fr),
    "info"                                => infostring,
  )
end

export rings_to_dict

function rings_to_dict(frs::Vector{FusionRing})
  infostring =
    """
    Dataset of fusion rings containing the fields "data", "info", and "order" with the following content.
    * data: contains JSON dictionary of fusion rings where the keys are UUID1 identifiers.
    * order: contains a list of identifiers of the fusion rings in the original order of the list of rings that was exported.
    * info: contains all info on how to interpret the data stored in this JSON file.

    Each fusion ring form the "data" is given by JSON dictionary.
    """ * ringfieldinfo

  frstrings = Dict(fr => fusion_ring_string(fr) for fr in frs)
  dt = []
  @showprogress dt = 0.5 desc="Converting rings to json object" for fr in frs
    push!(dt, frstrings[fr] => delete!(ring_to_dict(fr), "info"))
  end

  return Dict(
    "data"  => Dict(dt...),
    "info"  => infostring,
    "order" => [frstrings[fr] for fr in frs], # to preserve the order when importing
  )
end

function fusion_ring_string(fr::FusionRing)
  c = uuid(fr)
  if ismissing(c)
    id = string(UUIDs.uuid1())
    return id
  end
  return c
end

function fusion_ring_file_name(fr::FusionRing)
  awid = anyonwiki_code(fr)
  if !ismissing(awid)
    return join(string.(awid), "_") * ".json"
  else
    awid = uuid(fr)
  end

  !ismissing(awid) && return awid * ".json"

  return string(UUIDs.uuid1()) * ".json"
end

export export_ring

function export_ring(filename::String, fr::FusionRing)
  return write_json(filename, ring_to_dict(fr))
end

export export_rings

"""export_rings(fn::String,frs::Vector{FusionRing})

Exports the list of fusion rings frs. The list can be imported using import_rings
"""

function export_rings(filename::String, frs::Vector{FusionRing})
  d = rings_to_dict(frs)
  println("exporting json file")
  return write_json(filename, d)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                      ADDING NEW RINGS TO THE DATABASE                           ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export auto_complete_missing_info

"""auto_complete_missing_info(fr::FusionRing)::FusionRing

returns a fusion ring where the following data is computed and added if it was missing:
* characters
* sub_fusion_rings
* frobenius_perron_dimension(s)
* formal_codegrees
* all_gradings
* upper_central_series
* realizations (as tensor_product)
* automorphism_group
"""

function auto_complete_missing_info(fr::FusionRing)::FusionRing
  # is always defined
  mt = multiplication_table(fr)

  # should always be defined but just in case...
  id = ismissing(fr.uuid) ? UUIDs.uuid1() : fr.uuid

  # can't be autocompleted
  awc = fr.anyonwiki_code

  # autocompletion is too expensive on case-by-case basis
  nms = fr.names

  # autocompletion is too expensive on case-by-case basis
  tnms = fr.texnames

  # no autocompletion necessary
  lbls = fr.labels

  # automatically loads if stored and computes if missing
  chrs = is_commutative(fr) ? characters(fr) : missing

  sfr = sub_fusion_rings(fr)

  fpd = frobenius_perron_dimension(fr)

  fpds = frobenius_perron_dimensions(fr)

  fcds = formal_codegrees(fr)

  # way to expensive to autocomplete in general
  # TODO: can be (partially) completed for special cases
  hcwp = fr.has_categories_with_props

  # way to expensive to autocomplete in general
  cts = fr.categorifications

  # can't be autocompleted
  # TODO: if this software is used for certain properties the (future) corresponding paper should be added
  rfs = fr.references

  # can't be autocompleted
  # TODO: if this software is used for certain properties, it should be added
  sftw = fr.software

  # automatically loads if stored and computes if missing
  ags = all_gradings(fr)

  ucs = upper_central_series(fr)

  # need to manually check and add for this property
  rlztns = fr.realizations
  if ismissing(rlztns["tensor_product"]) || isnothing(rlztns["tensor_product"])
    rlztns["tensor_product"] = decompositions(fr; kind = "tensor_product")
  end

  # automatically loads if stored and computes if missing
  ag = automorphism_group(fr)

  return fusion_ring(
    mt;
    uuid                        = id,
    anyonwiki_code              = awc,
    names                       = nms,
    texnames                    = tnms,
    labels                      = lbls,
    characters                  = chrs,
    sub_fusion_rings            = sfr,
    frobenius_perron_dimension  = fpd,
    frobenius_perron_dimensions = fpds,
    formal_codegrees            = fcds,
    has_categories_with_props   = hcwp,
    categorifications           = cts,
    references                  = rfs,
    software                    = sftw,
    all_gradings                = ags,
    upper_central_series        = ucs,
    realizations                = rlztns,
    automorphism_group          = ag,
    skip_check                  = true, #Otherwise the input wouldn't be a FusionRing
  )
end

function auto_complete_missing_info(a::Array{Int64, 3})::FusionRing
  return auto_complete_missing_info(fusion_ring(a))
end

export change_properties

"""change_properties(fr::FusionRing,d::Dict{Symbol,Any})

takes a fusion ring fr and a dictionary d that maps field names of a FusionRing to new values and returns a new fusion ring where each field has that new value

change_properties(fr::FusionRing,:s => v)

is equal to change_property( fr, Dict( :s => v ) )
"""

function change_properties(fr::FusionRing, d::Dict{Symbol, T}) where {T}
  isempty(d) && return fr

  fns = fieldnames(FusionRing)

  kys = collect(keys(d))

  invalid_keys = setdiff(kys, fns)
  !isempty(invalid_keys) &&
    error("The following keys are invalid fields of a FusionRing: $invalid_keys")

  vals = []

  for fn in fns
    if fn ∈ kys
      push!(vals, d[fn])
    else
      push!(vals, getproperty(fr, fn))
    end
  end

  return FusionRing(vals...)
  #result = fr
  #for (k,v) in pairs(d)
  #  @info (k,v)
  #  result = @set result.$k = v
  #end
  #return result
  #
end

function change_properties(fr::FusionRing, tup::Pair{Symbol, T}) where {T}
  return change_properties(fr, Dict(tup))
end

export add_names

function add_names(fr::FusionRing, tup::Pair{String, Vector{String}})
  return add_nms(fr, tup; kind = "nontex")
end

export add_tex_names

function add_tex_names(fr::FusionRing, tup::Pair{String, Vector{String}})
  return add_nms(fr, tup; kind = "tex")
end

function add_nms(fr::FusionRing, tup::Pair{String, Vector{String}}; kind = "tex")
  id, newnames = tup
  if id ∉ _naming_priority_order
    error("$(tup[1]) should be an element of $_naming_priority_order")
  end

  nms = kind == "tex" ? deepcopy(tex_names(fr)) : deepcopy(names(fr))

  if ismissing(nms[id])
    nms[id] = newnames
  else
    nms[id] = vcat(newnames, nms[id])
  end

  if kind == "tex"
    return change_properties(fr, :texnames => nms)
  else
    return change_properties(fr, :names => nms)
  end
end
