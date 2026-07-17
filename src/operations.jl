
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
  return sort( [c for (c, m) in fusion_product(fr, a, b) if m>0] )
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
    if p ≠ nothing
      target = ring
      break
    end
  end

  if p === nothing
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

  # TODO: might want to conserve as much information as possible, but 
  # it's better to wait until all fields of the FusionRing struct are finalized
  return fusion_ring(Nsub)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    permute                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: implement permute for permutations from OSCAR
#fix: permute changes meaning of every simple-object index
#previously: multp. table was permuted but data attached to indices was left same
#mental model for recollection (for me lol):
#ring: labels = ["0", "1", "2"]
#perm = [1, 3, 2]
#new_mt[a, b, c] =
# old_mt[perm[a], perm[b], perm[c]]
"""
    permute(perm, ring) -> FusionRing

Return  isomorphic presentation of `ring` in which new simple object `i`
corresponds to old simple object `perm[i]`.

 permutation must fix the tensor unit.
"""
function permute(
  perm::Vector{Int},
  ring::FusionRing,
)::FusionRing
  r = rank(ring)

  length(perm) == r ||
    throw(ArgumentError("permutation length must equal the ring rank"))

  sort(perm) == collect(1:r) ||
    throw(ArgumentError("perm must contain every index from 1 to $r exactly once"))

  perm[1] == 1 ||
    throw(ArgumentError("the tensor unit at index 1 must remain fixed"))

  new_mt = permute_mult_tab(
    multiplication_table(ring),
    perm,
  )

  result = @set ring.multiplication_table = new_mt

  # New id i -> old index perm[i].
  result = @set result.labels = ring.labels[perm]

  # Exact FP dims are indexed by simple objects.
  if ring.frobenius_perron_dimensions !== missing
    result = @set result.frobenius_perron_dimensions =
      ring.frobenius_perron_dimensions[perm]
  end

  # Numeric FP dims  also indexed by simple objects.
  if ring.numeric_frobenius_perron_dimensions !== missing
    result = @set result.numeric_frobenius_perron_dimensions =
      ring.numeric_frobenius_perron_dimensions[perm]
  end

  # char table cols correspond to simple objects.
  if ring.characters !== missing
    result = @set result.characters =
      ring.characters[:, perm]
  end

  if ring.numeric_characters !== missing
    result = @set result.numeric_characters =
      ring.numeric_characters[:, perm]
  end

  #  every old subring index  its new position.
  if ring.sub_fusion_rings !== missing
    old_to_new = invperm(perm)

    permuted_subrings = [
      begin
        new_subring = copy(subring)
        new_subring["injection"] =
          old_to_new[subring["injection"]]
        new_subring
      end
      for subring in ring.sub_fusion_rings
    ]

    result = @set result.sub_fusion_rings =
      permuted_subrings
  end

  #  representations are indexed by simple objects, but  precise
  # stored schemas need separate transformation methods. Keeping them would
  # attach incorrectly indexed data to the permuted ring.
  if ring.projective_SL2Z_reps !== missing
    result = @set result.projective_SL2Z_reps = missing
  end

  if ring.numeric_projective_SL2Z_reps !== missing
    result = @set result.numeric_projective_SL2Z_reps = missing
  end

  return result
end

"""
    permute_mult_tab(N, p)

Apply permutation `p` (fixing 1) to all three indices of `N`.
"""
#fix: replaced triple loop with direct indexing
function permute_mult_tab(
  table::Array{Int,3},
  perm::Vector{Int},
)::Array{Int,3}
  r = size(table, 1)

  size(table) == (r, r, r) ||
    throw(ArgumentError("multiplication table must be cubic"))

  length(perm) == r ||
    throw(ArgumentError("permutation length must equal table rank"))

  sort(perm) == collect(1:r) ||
    throw(ArgumentError("perm is not a valid permutation"))

  perm[1] == 1 ||
    throw(ArgumentError("permutation must fix the tensor unit"))

  return table[perm, perm, perm]
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
  sdlist = reduce( vcat, sd, init = Int[] )[2:end]
  sort!(sdlist; by = i -> fpd[i], rev = (order == :decreasing))

  # sort pairs of dual elements and flatten the list
  sort!(nsd; by = p -> fpd[p[1]], rev = (order == :decreasing))
  nsdlist = reduce(vcat, nsd; init = Int[])

  return vcat( [1], sdlist, nsdlist)
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

  names_tp = if (isempty(names(r1)) || isempty(names(r2)))
    String[]
  else
    [string(names(r1)[1], "⊗", names(r2)[1])]
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

function tensor_product(rings::Vector{FusionRing})::FusionRing
  isempty(rings) && error("Need at least one fusion ring")
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
