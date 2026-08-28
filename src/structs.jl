export FusionRing

struct FusionRing
  multiplication_table::Array{Int, 3}

  uuid::Union{Missing, String}

  anyonwiki_code::Union{Missing, Vector{Int64}}

  names::Dict{String, Union{Missing, Vector{String}}}

  texnames::Dict{String, Union{Missing, Vector{String}}}

  labels::Array{String, 1}

  characters::Union{Missing, Matrix{String}}

  sub_fusion_rings::Union{Missing, Vector{Tuple{Vector{Int64}, String}}}

  frobenius_perron_dimension::Union{Missing, String}

  frobenius_perron_dimensions::Union{Missing, Vector{String}}

  formal_codegrees::Union{Missing, Vector{String}}

  has_categories_with_props::Dict{String, Vector{Any}}

  categorifications::Union{Missing, Vector{String}}

  references::Union{Missing, Dict{String, Vector{String}}}

  software::Union{Missing, Dict{String, Vector{String}}}

  all_gradings::Union{Missing, Vector{Tuple{Vector{Int64}, String}}}

  upper_central_series::Union{Missing, Vector{Tuple{Vector{Int64}, String}}}

  realizations::Union{Missing, Dict{String, Any}}

  automorphism_group::Union{Missing, PermGroup}
end

export fusion_ring

function check_struct_const(mt)
  return all(x -> x isa Integer && x >= 0, mt)
end

function check_mt_dims(mt)
  dims = size(mt)
  return length(dims) == 3 && first(dims) >= 1 && is_constant_array(dims)
end

function check_unit(mt)
  r = size(mt)[1]
  δ(i, j) = i == j ? 1 : 0
  for i in 1:r, j in 1:r
    if !(mt[1, i, j] == mt[i, 1, j] == δ(i, j))
      return false
    end
    continue
  end
  return true
end

function check_inverse(mt)
  r = size(mt, 1)
  unit_coefficients = @view mt[:, :, 1]

  for i in 1:r
    left_coefficients = view(unit_coefficients, :, i)
    right_coefficients = view(unit_coefficients, i, :)
    left_duals = findall(==(1), left_coefficients)
    right_duals = findall(==(1), right_coefficients)

    length(left_duals) == 1 || return false
    length(right_duals) == 1 || return false
    only(left_duals) == only(right_duals) || return false

    all(x -> x == 0 || x == 1, left_coefficients) || return false
    all(x -> x == 0 || x == 1, right_coefficients) || return false
  end

  return true
end

function check_associativity(mt::Array{Int, 3})
  r = size(mt, 1)
  for a in 1:r, b in 1:r, c in 1:r, d in 1:r
    lhs = sum(mt[a, b, e] * mt[e, c, d] for e in 1:r)
    rhs = sum(mt[a, f, d] * mt[b, c, f] for f in 1:r)
    lhs == rhs || return false
  end
  return true
end

function check_labels(mt, names)
  return length(names) == size(mt, 1)
end

function is_valid_uuid(s::String)
  return !isnothing(tryparse(UUID, s))
end

function fusion_ring(
  mt::Array{Int, 3};
  uuid                        = missing,
  anyonwiki_code              = missing,
  names                       = missing,
  texnames                    = missing,
  labels                      = missing,
  characters                  = missing,
  sub_fusion_rings            = missing,
  frobenius_perron_dimension  = missing,
  frobenius_perron_dimensions = missing,
  formal_codegrees            = missing,
  has_categories_with_props   = missing,
  categorifications           = missing,
  references                  = missing,
  software                    = missing,
  all_gradings                = missing,
  upper_central_series        = missing,
  realizations                = missing,
  automorphism_group          = missing,
  skip_check                  = false,
)
  if !skip_check
    check_struct_const(mt) || error("All structure constants must be non-negative integers")
    check_mt_dims(mt) ||
      error("multiplication_table must be a 3-tensor with equal side lengths")
    check_unit(mt) || error("First basis element must act as unit object")
    check_inverse(mt) || error("Each simple object must have a unique inverse")
    check_associativity(mt) || error("Structure constants violate associativity")
    if labels !== [] && !ismissing(labels) && !check_labels(mt, labels)
      error("labels length ≠ rank")
    end
  end

  #@info (mt::Array{Int,3},  uuid, anyonwiki_code, names, texnames, labels, characters, sub_fusion_rings, frobenius_perron_dimension, frobenius_perron_dimensions, formal_codegrees, has_categories_with_props           ,  categorifications                   ,  references                          ,  software                            ,  all_gradings                        ,  upper_central_series,  realizations, automorphism_group)

  (ismissing(labels) || labels == []) &&
    (labels = String[bold_integer(i) for i in 1:size(mt, 1)])

  id = ismissing(uuid) ? string(UUIDs.uuid1()) : uuid

  !is_valid_uuid(id) && error("UUID string is invalid")

  r   = size(mt, 1)
  sd  = count(x -> x == 1, diag(mt[:, :, 1]))
  nsd = r - sd

  fpd = if frobenius_perron_dimension isa Oscar.QQBarFieldElem
    qqb_id(frobenius_perron_dimension)
  else
    frobenius_perron_dimension
  end

  fpds = if frobenius_perron_dimensions isa Vector{Oscar.QQBarFieldElem}
    qqb_id.(frobenius_perron_dimensions)
  else
    frobenius_perron_dimensions
  end

  fcds = if formal_codegrees isa Vector{Oscar.QQBarFieldElem}
    qqb_id.(formal_codegrees)
  else
    formal_codegrees
  end

  # characters can be encoded in a variety of different ways
  mtspceqqb = AbstractAlgebra.Generic.MatSpaceElem{Nemo.QQBarFieldElem}

  chars = if characters isa mtspceqqb || characters isa Matrix{QQBarFieldElem}
    [qqb_id(characters[i, j]) for i in 1:r, j in 1:r]
  elseif characters isa Vector{Vector{String}}
    [characters[i][j] for i in 1:r, j in 1:r]
  elseif characters isa Vector{Vector{Any}}
    String[characters[i][j] for i in 1:r, j in 1:r]
  else
    characters
  end

  function convert_nms(x)
    ismissing(x) && return missing
    isnothing(x) && return missing
    x isa Vector{String} && return x
    all(el -> el isa String, x) && return string.(x)

    return error(
      "Values of names dictionary should be nothing, missing or a list of Strings"
    )
  end

  nms = if ismissing(names)
    _nonames
  elseif names isa Vector
    mscnames(string.(names))
  else
    Dict(k => convert_nms(v) for (k, v) in names)
  end

  texnms = if ismissing(texnames)
    _nonames
  elseif texnames isa Vector
    mscnames(string.(names))
  else
    Dict(k => convert_nms(v) for (k, v) in texnames)
  end

  # LEGACY compatibility
  hcwp = if has_categories_with_props isa Vector
    vec = has_categories_with_props
    Dict{String, Vector{Any}}(row[1] => [row[2], row[3][1], row[3][2]] for row in vec)
  elseif ismissing(has_categories_with_props)
    Dict(
      "Fusion"    => [missing, "", "Unknown to AnyonWiki"],
      "Pivotal"   => [missing, "", "Unknown to AnyonWiki"],
      "Unitary"   => [missing, "", "Unknown to AnyonWiki"],
      "Spherical" => [missing, "", "Unknown to AnyonWiki"],
      "Braided"   => [missing, "", "Unknown to AnyonWiki"],
      "Ribbon"    => [missing, "", "Unknown to AnyonWiki"],
      "Modular"   => [missing, "", "Unknown to AnyonWiki"],
    )
  else
    has_categories_with_props
  end

  #TODO: let references be paper of FusionRings if missing.
  refs = if references isa Vector
    Dict{String, Vector{String}}("All" => references)
  else
    references
  end

  #TODO: let software be current version of FusionRings if missing.
  sftw = if software isa Vector
    Dict{String, Vector{String}}("All" => software)
  else
    software
  end

  ag = if all_gradings isa Vector{Tuple{Vector{Int64}, FusionRing}}
    [(t[1], t[2].uuid) for t in all_gradings]
  else
    all_gradings
  end

  ucs = if upper_central_series isa Vector{Tuple{Vector{Int64}, FusionRing}}
    [(t[1], t[2].uuid) for t in upper_central_series]
  else
    upper_central_series
  end

  sfr = if sub_fusion_rings isa Vector{Tuple{Vector{Int64}, FusionRing}}
    [(t[1], t[2].uuid) for t in sub_fusion_rings]
  else
    sub_fusion_rings
  end

  rlztns =
    ismissing(realizations) ? Dict{String, Any}("tensor_product" => missing) : realizations

  return FusionRing(
    mt,
    id,
    anyonwiki_code,
    nms,
    texnms,
    labels,
    chars,
    sfr,
    fpd,
    fpds,
    fcds,
    hcwp,
    categorifications,
    refs,
    sftw,
    ag,
    ucs,
    rlztns,
    automorphism_group,
  )
end
