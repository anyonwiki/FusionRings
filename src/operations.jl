
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                          various helper functions                               ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""
    fusion_product(fr, a, b) -> Dict{Int,Int}

Return the decomposition of `a ⊗ b` as a multiplicity dictionary
`Dict{simple_index => multiplicity}`.
"""
function fusion_product(fr::FusionRing, a::Int, b::Int)
  N = @views multiplication_table(fr)[a, b, :]
  out = Dict{Int, Int}()
  @inbounds for (c, m) in enumerate(N)
    m==0 && continue
    out[c] = m
  end
  return out
end

"Return vector of simple indices with positive multiplicity in `a × b`."
function fusion_outcomes(fr::FusionRing, a::Int, b::Int)::Vector{Int}
  return sort([c for (c, m) in fusion_product(fr, a, b) if m>0])
end

# TODO: implement change_property function using the Accessors package

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                replace_by_known                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export replace_by_known

# keep_order=true permutes the elements of the known ring to match the
# elements of the given ring.
# safe_return=true makes the function return the original ring in case no built-in 
# ring is found

function replace_by_known(fr::FusionRing; keep_order = true, safe_return = true)::FusionRing
  r = rank(fr)
  m = multiplicity(fr)
  n = nnsd(fr)

  if !haskey(frd, [r, m, n])
    safe_return && return fr
    return missing
  end

  proposals = collect(values(frd[[r, m, n]]))

  # find the permutation
  p      = nothing
  target = nothing
  for ring in proposals
    p = which_permutation(ring, fr)
    if !isnothing(p)
      target = ring
      break
    end
  end

  if isnothing(p)
    safe_return && return fr
    return missing
  end

  if keep_order
    permute(p[1], target)
  else
    return target
  end
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                restrict_subring                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export restrict_subring

# restrict ring to subindices S (Anyonica's MT[ring][[el,el,el]])
function restrict_subring(
  fr::FusionRing, S::Vector{Int}; check_closed::Bool = true
)::FusionRing
  sS = sort(S)

  sS == collect(1:rank(fr)) && return fr

  unique(sS) != sS && error("vector should contain each element only once")

  N = multiplication_table(fr)
  Nsub = N[sS, sS, sS]

  if check_closed
    # sanity: S must be fusion-closed
    rmap = zeros(Int, rank(fr))
    @inbounds for (k, v) in enumerate(sS)
      rmap[v] = k
    end
    @inbounds for a in sS, b in S
      for c in findall(>(0), N[a, b, :])
        rmap[c] != 0 || error("subset not fusion-closed")
      end
    end
  end

  return fusion_ring(Nsub; labels = labels(fr)[sS])
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    permute                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: implement permute for permutations from OSCAR
export permute

"""permute(perm::Vector{Int},r::FusionRing)::FusionRing – return a new `FusionRing` with all data
    permuted by `perm`.
   `perm[1]` **must** equal 1 to keep the vacuum first.
"""

function permute(σ::Vector{Int}, r::FusionRing)::FusionRing
  n = rank(r)
  n == length(σ) || throw(ArgumentError("perm length ≠ rank"))
  sort(σ) == collect(1:n) || throw(ArgumentError("perm must be a true permutation"))
  σ[1] == 1 || throw(ArgumentError("vacuum must stay at index 1"))
  σ == collect(1:n) && return r

  pmt = permute_mult_tab(multiplication_table(r), σ)
  iσ = invperm(σ)
  permute_indices(v) = iσ[v]

  lbs = r.labels[σ]

  # Character columns and FP dimensions are indexed by simple objects.
  ch = r.characters
  if !ismissing(ch)
    ch = ch[:, σ]
  end

  sr = r.sub_fusion_rings
  if !ismissing(sr)
    sr = [(permute_indices(t[1]), t[2]) for t in sr]
  end

  fpd = r.frobenius_perron_dimensions
  if !ismissing(fpd)
    fpd = fpd[σ]
  end

  ag = r.all_gradings
  if !ismissing(ag)
    ag = [(g[1][σ], g[2]) for g in ag]
  end

  aut = r.automorphism_group
  if !ismissing(aut)
    aut = conjugate_group(aut, perm(iσ))
  end

  ucs = r.upper_central_series
  if !ismissing(ucs)
    ucs = [(permute_indices(t[1]), t[2]) for t in ucs]
  end

  return change_properties(
    r,
    Dict{Symbol, Any}(
      :multiplication_table        => pmt,
      :labels                      => lbs,
      :characters                  => ch,
      :sub_fusion_rings            => sr,
      :frobenius_perron_dimensions => fpd,
      :upper_central_series        => ucs,
      :all_gradings                => ag,
      :automorphism_group          => aut,
    ),
  )
end

"""
    permute_mult_tab(N, p)

Apply permutation `p` (fixing 1) to all three indices of `N`.
"""
function permute_mult_tab(N::Array{Int, 3}, p::Vector{Int})::Array{Int, 3}
  r = size(N, 1)
  size(N) == (r, r, r) || throw(ArgumentError("multiplication table must be cubic"))
  length(p) == r || throw(ArgumentError("permutation length must equal table rank"))
  sort(p) == collect(1:r) || throw(ArgumentError("perm must be a true permutation"))
  p[1] == 1 || throw(ArgumentError("permutation must fix the unit at index 1"))

  return N[p, p, p]
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                     sort                                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: figure out which one is best: string or symbol as named input 
#TODO: make sort function have secondary sort method based on e.g. grouping dual 
# elements together

export sort

function sort(fr::FusionRing; by = "fpdims", order::Symbol = :increasing)
  if by == "fpdims"
    return permute(perm_vec_qd(fr; order = order), fr)
  elseif by == "sd_conj"
    return permute(perm_vec_sd_conj(fr; order = order), fr)
  else
    message("by= argument was not \"fpdims\" or \"sd_conj\".")
  end
end

"""perm_vec_qd(r; order = :increasing) – permutation that sorts the non‑vacuum
    elements by Frobenius–Perron dimension."""
function perm_vec_qd(r::FusionRing; order::Symbol = :increasing)::Vector{Int}
  idx = collect(2:rank(r))
  qd  = fpdims(r)
  sort!(idx; by = i -> qd[i], rev = (order == :decreasing))
  return vcat(1, idx)
end

"""perm_vec_sd_conj(r; order = :increasing) – self‑duals first (ordered by fpdim), then conjugate
    pairs, with blocks ordered by FP‑dimension."""
function perm_vec_sd_conj(r::FusionRing; order::Symbol = :increasing)::Vector{Int}
  fpd   = fpdims(r)
  pairs = conjugate_pairs(r)

  sd, nsd = binsplit(p -> length(p)==1, pairs)

  # flatten list of self-dual elements, remove unit, and and sort by fpdim
  sdlist = reduce(vcat, sd; init = Int[])[2:end]
  sort!(sdlist; by = i -> fpd[i], rev = (order == :decreasing))

  # sort pairs of dual elements and flatten the list
  sort!(nsd; by = p -> fpd[p[1]], rev = (order == :decreasing))
  nsdlist = reduce(vcat, nsd; init = Int[])

  return vcat([1], sdlist, nsdlist)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                tensor_product                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export tensor_product

function tensor_product(r1::FusionRing, r2::FusionRing)::FusionRing
  m, n = rank(r1), rank(r2)
  mt1, mt2 = multiplication_table(r1), multiplication_table(r2)

  mt = zeros(Int, m*n, m*n, m*n)
  @inbounds for a in 1:m, α in 1:n, b in 1:m, β in 1:n, c in 1:m, γ in 1:n
    i = (a-1)*n + α
    j = (b-1)*n + β
    k = (c-1)*n + γ
    mt[i, j, k] = mt1[a, b, c] * mt2[α, β, γ]
  end

  # Assemble element names
  elnames = [string(e1, "⊗", e2) for e1 in labels(r1) for e2 in labels(r2)]

  names_tp = if ismissing(name(r1)) || ismissing(name(r2))
    _nonames
  else
    tpnames(names(r1), names(r2))
  end

  fpdims_new = try
    [d1 * d2 for d1 in fpdims(r1) for d2 in fpdims(r2)]
  catch
    missing
  end

  return fusion_ring(
    mt; names = names_tp, labels = elnames, frobenius_perron_dimensions = fpdims_new
  )
end

function tpnames(nms1, nms2)
  d = Dict{String, Union{Missing, Vector{String}}}()
  keys1 = keys(nms1)
  keys2 = keys(nms2)
  sort(collect(keys1)) != sort(collect(keys2)) &&
    error("The naming keys of the fusion rings do not match")

  for k in keys1
    nk1 = nms1[k]
    nk2 = nms2[k]
    if !ismissing(nk1) && !ismissing(nk2)
      d[k] = [n1 * "⊗" * n2 for n1 in nk1 for n2 in nk2]
    else
      d[k] = missing
    end
  end

  # add products of names of different conventions to miscellaneous

  for k1 in keys1, k2 in keys2
    # we already covered k1 === k2
    k1 === k2 && continue

    nk1 = nms1[k1]
    nk2 = nms2[k2]
    im1 = ismissing(nk1)
    im2 = ismissing(nk2)
    if !im1 && !im2
      if ismissing(d["miscellaneous"])
        d["miscellaneous"] = String[]
      end

      for n1 in nk1, n2 in nk2
        push!(d["miscellaneous"], n1 * "⊗" * n2)
      end
    end
  end

  return d
end

function tensor_product(rings::Vector{FusionRing})::FusionRing
  # empty product is unit by default
  isempty(rings) && return from_anyonwiki_code([1, 1, 0, 1])
  length(rings) == 1 && return rings[1]

  out = rings[1]
  for R in rings[2:end]
    out = tensor_product(out, R)
  end
  return out
end

function tensor_product(rings...)::FusionRing
  return tensor_product([rings...])
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   to_group                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export to_group

function to_group(fr::FusionRing)
  ct = cayley_table(fr)
  gm(a, b) = ct[a, b]
  return permutation_group(generic_group(1:rank(fr), gm)[1])
end

# TODO: implement bicrossed product. @Szagha02: not a priority 
# @gvercley will do this at some point
