
#function change_fusion_ring_property(r::FusionRing, dict)

#end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            multiplication_table                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export multiplication_table

function multiplication_table(r::FusionRing)::Array{Int, 3}
  return r.multiplication_table
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    rank                                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: write test: rank(ring) = anyonwiki_code(ring)[1]

export rank

function rank(r::FusionRing)::Int
  return size(multiplication_table(r), 1)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    names                                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export names

function names(r::FusionRing)::Array{String, 1}
  return r.names
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  tex_names                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export tex_names

function tex_names(r::FusionRing)::Array{String, 1}
  return r.texnames
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    labels                                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export labels

function labels(r::FusionRing)::Array{String, 1}
  return r.labels
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                               conjugation_matrix                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export conjugation_matrix

function conjugation_matrix(fr::FusionRing)
  @views multiplication_table(fr)[:, :, 1]
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  multiplicity                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#TODO: write test: multiplicity(r) = anyonwiki_code(r)[2]

export multiplicity

function multiplicity(r::FusionRing)::Int
  return maximum(multiplication_table(r))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      mult                                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export mult

mult = multiplicity

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           nonzero_structure_constants                           ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: write test: check whether number of nonzero struct const is same as for data in
# folder test/testdata/properties/

export nonzero_structure_constants

function nonzero_structure_constants(r::FusionRing)::Vector{Tuple{Int64, Int64, Int64}}
  mt = multiplication_table(r)
  return Tuple.(findall(x -> x > 0, mt))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      nzsc                                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export nzsc

nzsc = nonzero_structure_constants

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           frobenius_perron_dimensions                           ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export frobenius_perron_dimensions

function frobenius_perron_dimensions(
  r::FusionRing; force_compute = false
)::Vector{QQBarFieldElem}
  stored_dims = r.frobenius_perron_dimensions
  if stored_dims===missing || force_compute
    mt = multiplication_table(r)
    multmats = [matrix(ZZ, mt[i, :, :]) for i in 1:rank(r)]
    return [first(eigenvalues(QQBar, A)) for A in multmats]
  else
    return stored_dims
  end
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                     fpdims                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export fpdims

fpdims = frobenius_perron_dimensions

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           frobenius_perron_dimension                            ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export frobenius_perron_dimension

function frobenius_perron_dimension(r::FusionRing)::QQBarFieldElem
  return sum(fpdims(r) .^ 2)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                     fpdim                                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export fpdim

fpdim = frobenius_perron_dimension

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  numeric_fpdims                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export numeric_fpdims

function numeric_fpdims(fr::FusionRing)
  r = rank(fr)
  S = zeros(Float64, r, r)
  N = multiplication_table(fr)
  for a in 1:r
    @views S .+= N[a, :, :]
  end
  vals, vecs = eigen(S)
  idx = argmax(vals)
  v = abs.(vecs[:, idx])
  return v ./ v[1]
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   numeric_fpdim                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export numeric_fpdim

function numeric_fpdim(fr::FusionRing)
  return sum(x->x*x, numeric_fpdims(fr))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           num_self_dual_non_self_dual                           ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export num_self_dual_non_self_dual

function num_self_dual_non_self_dual(r::FusionRing)::Array{Int, 1}
  sd  = count(x -> x == 1, diag(conjugation_matrix(r)))
  nsd = rank(r) - sd
  return Int[sd, nsd]
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                     nsdnsd                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export nsdnsd

nsdnsd = num_self_dual_non_self_dual

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  num_self_dual                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export num_self_dual

function num_self_dual(r::FusionRing)::Int
  return first(nsdnsd(r))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      nsd                                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export nsd

nsd = num_self_dual

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                num_non_self_dual                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: write test: num_non_self_dual(ring) == anyonwiki_code(ring)[3]
export num_non_self_dual

function num_non_self_dual(r::FusionRing)::Int
  return last(nsdnsd(r))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                     nnsd                                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export nnsd

nnsd = num_non_self_dual

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  is_group_ring                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export is_group_ring

function is_group_ring(r::FusionRing)::Bool
  return sum(multiplication_table(r)) == rank(r)^2
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   cayley_table                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export cayley_table

function cayley_table(fr::FusionRing)
  !is_group_ring(fr) && message("Ring must be group ring")

  mt = multiplication_table(fr)
  r  = rank(fr)
  return [findfirst(==(1), mt[a, b, :]) for a in 1:r, b in 1:r]
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                conjugate_element                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export conjugate_element

function conjugate_element(r::FusionRing)
  return a -> conjugate_element(r, a)
end

"""
    conjugate_element(fr, a) -> Int

Return the integer index of the dual (conjugate) simple object of `a`.
Accepts an integer index, a `String`, or a `Symbol`.
"""
function conjugate_element(fr::FusionRing, a::Int64)::Int64
  C = conjugation_matrix(fr)
  return findfirst(==(1), C[a, :])
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                 conjugate_pairs                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export conjugate_pairs

function conjugate_pairs(fr::FusionRing)::Vector{Vector{Int}}
  d  = conjugate_element(fr);
  us = unique ∘ sort
  return unique([us([a, d(a)]) for a in 1:rank(fr)])
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                 anyonwiki_code                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export anyonwiki_code

function anyonwiki_code(r::FusionRing)::Union{Array{Int, 1}, Missing}
  return r.anyonwiki_code
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    barcode                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export barcode

function barcode(r::FusionRing)
  return r.barcode
end

function mult_tab_code(mat::Array{Int, 2}, mult::Int)::Int
  return error("mult_tab_code not implemented yet")
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                               sub_fusion_rings                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export sub_fusion_rings

function sub_fusion_rings(r::FusionRing)
  dictvec = r.sub_fusion_rings
  if dictvec !== missing
    [
      Dict("injection"   => dict["injection"], "fusion_ring" => awc(dict["anyonwiki_code"]))
      for dict in dictvec
    ]
  else
    subsets = sub_fusion_ring_subsets
    [
      Dict("injection"   => s, "fusion_ring" => replace_by_known(fusion_ring(mt[s, s, s])))
      for s in subsets
    ]
  end
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                             sub_fusion_ring_subsets                             ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export sub_fusion_ring_subsets

"""
    sub_fusion_ring_subsets(fr::FusionRing) -> Vector{Vector{Int}}

Enumerate all **proper, nontrivial** fusion-closed subsets of elements containing
the unit (index 1), returned as **index vectors**.
"""

function sub_fusion_ring_subsets(fr::FusionRing)::Vector{Vector{Int}}
  r = rank(fr)
  r <= 2 && return Vector{Int}[]

  result = Vector{Int}[]

  # group elements by conjugacy since if a is part of subring a* also has to be
  pairs = conjugate_pairs(fr)[2:end]
  np    = size(pairs, 1)

  # flatten conjugate_pairs back to 1D vector
  function flatten(l::Vector{Vector{Int}})
    res = Int[]
    for vec in l, el in vec
      push!(res, el)
    end
    return res
  end

  for k in 1:(np - 1), subset in combinations(pairs, k)
    S = vcat([1], flatten(collect(subset)))
    if is_sub_fusion_ring(fr, S)
      push!(result, S)
    end
  end
  return result
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              is_sub_fusion_ring                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export is_sub_fusion_ring

#TODO: write test that checks whether all fusion rings obained via sub_fusion_rings are fusion rings
#TODO: write test: pick some rings that are not sub fusion rings and check whether this is recognized

"""
    is_sub_fusion_ring(fr, S) -> Bool

Return `true` iff `S` is a fusion-closed subset of simples containing the unit.

`S` may be a vector of indices (`Int`/`Integer`) or a vector of labels
(`String`/`Symbol`).
"""
function is_sub_fusion_ring(fr::FusionRing, S::Vector{Int})::Bool
  # any subring must contain unit
  isempty(S) && return false
  1 ∉ S && return false

  # indices in S must lie in range 1, ...,  r
  r = rank(fr)
  !all(i -> 1 <= i <= r, S) && return false

  mt = multiplication_table(fr)[S, S, S]

  return check_struct_const(mt) &&
           check_mt_dims(mt) &&
           check_unit(mt) &&
           check_inverse(mt) &&
           check_associativity(mt)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            is_equivalent_fusion_ring                            ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#TODO: write test: permute some rings and check whether is_equivalent_fusion_ring returns true
#TODO: write test: take some rings with different codes and check whether we get false

export is_equivalent_fusion_ring

function is_equivalent_fusion_ring(ring1::FusionRing, ring2::FusionRing)::Bool
  return which_permutation(ring1, ring2) !== nothing
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                which_permutation                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#TODO: write test: permute ring using a vector and check whether which_permutation(..., all=true)
# containts this vector

export which_permutation

"""
    which_permutation(fr1, fr2, all=false, short_circuit=true)::Union{Nothing,Vector{Int64}}

    Returns a list with a permutation such that `permute( perm, fr1 ) == fr2` or nothing if none exist.
    If all=true, returns all permutations such that `permute( perm, fr1 ) == fr2` or nothing if none exist.
    If short_circuit=false, doesn't compute invariants to determine whether rings are
    not equivalent. This saves time if it is known a priori that the rings have to be equivalent.
"""
function which_permutation(
  fr1::FusionRing, fr2::FusionRing; all = false, short_circuit = true
)::Union{Nothing, Vector{Vector{Int64}}}
  r = rank(fr1)

  # if ranks are different we can't compare so necessary to short_circuit
  # since computation of number selfdual & non-selfdual elements is as
  # fast and gives stronger invariant, we use that one instead
  nsdnsd(fr1) != nsdnsd(fr2) && return nothing

  r == 1 && return [[1]]

  # check whether rings have same multiplicity
  short_circuit && mult(fr1) != mult(fr2) && return nothing

  # rings must have same fpdims
  fpd1 = fpdims(fr1)
  fpd2 = fpdims(fr2)

  short_circuit && sort(fpd1) ≠ sort(fpd2) && return nothing

  # check whether number of fusion outcomes per multiplicity
  # is the same for each element of fr1 and fr2
  m1 = multiplication_table(fr1)
  m2 = multiplication_table(fr2)

  grp1 = diag_channel_count(m1)
  grp2 = diag_channel_count(m2)

  short_circuit && sort(grp1) ≠ sort(grp2) && return nothing

  # construct invariants for all elements of both rings
  # these are triples of their fusion outcome count, fpdims
  # and being self-conjugate or not
  # we remove the first element since it needs to be fixed

  sc1 = is_self_conjugate(fr1).(collect(1:r))
  sc2 = is_self_conjugate(fr2).(collect(1:r))

  inv1 = collect(zip(grp1, fpd1, sc1))[2:end]
  inv2 = collect(zip(grp2, fpd2, sc2))[2:end]

  # for some operations the unit needs to be added again
  addunit(perm_vec) = vcat(1, perm_vec .+ 1)

  #sort invA and invB such that equal elements are next to each other.
  σ1 = sortperm(inv1)
  σ2 = sortperm(inv2)

  #permute the mult tabs so their invariants are sorted
  uσ1 = addunit(σ1)
  uσ2 = addunit(σ2)
  sm1 = m1[uσ1, uσ1, uσ1]
  sm2 = m2[uσ2, uσ2, uσ2]

  sorted_inv1 = inv1[σ1]

  S, _ = symmetries(sorted_inv1; sorted = true)

  if !all  # only need first permutation
    for σ in S
      p = addunit(Vector(σ, r-1))
      if sm1[p, p, p] == sm2
        iuσ2 = invperm(uσ2)
        return [iuσ2[p[uσ1]]]
      end
    end
  else # want all permutations
    allperms = Vector{Int}[]

    for σ in S
      p = addunit(Vector(σ, r-1))
      if sm1[p, p, p] == sm2
        iuσ2 = invperm(uσ2)
        push!(allperms, iuσ2[p[uσ1]])
      end
    end

    return allperms
  end
end

# returns the symmetry group S of the vector v together with
# the permutation σv that sorts v. The symmetries are thus of the form
# inverseperm(σv) ∘ g ∘ σv. If sorted=true it is assumed that v is sorted

function symmetries(v; sorted = false)::Tuple{PermGroup, PermGroupElem}
  is_empty(v) && return nothing

  n = size(v, 1)

  tocycles(v) = perm(symmetric_group(n), v)

  n == 1 && return (symmetric_group(1), tocycles([1]))

  if sorted
    return (_sorted_symmetries(v), tocycles(1:n))
  else
    σv = sortperm(v)
    vs = v[σv]
    return (_sorted_symmetries(vs), tocycles(σv))
  end
end

function _sorted_symmetries(v)::PermGroup
  return inner_direct_product(symmetric_group.(tally(v)[2]))
end

function diag_channel_count(N::Array{Int, 3})::Vector{Tuple{Vector{Int}, Vector{Int}}}
  return [tally(N[i, i, :]; sort = true) for i in 1:size(N, 1)]
end

function is_self_conjugate(fr)
  return x -> (x == conjugate_element(fr, x))
end

# Apply permutation P on all three indices: A'[i,j,k] = A[P[i],P[j],P[k]]
function _permute_multtab(A::Array{Int, 3}, P::Vector{Int})::Array{Int, 3}
  r = size(A, 1)
  B = similar(A)
  @inbounds for i in 1:r, j in 1:r, k in 1:r
    B[i, j, k] = A[P[i], P[j], P[k]]
  end
  return B
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            fusion_ring_automorphisms                            ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export fusion_ring_automorphisms
"""
    fusion_ring_automorphisms(fr) -> Vector{Vector{Int}}

Return all permutations `p` whose action on the indices leave the structure constants invariant.
"""
#TODO: write test: should match data in test/testdata/properties

function fusion_ring_automorphisms(fr::FusionRing)
  return which_permutation(fr, fr; all = true)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  decompositions                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export decompositions

function decompositions(fr::FusionRing, product = "TensorProduct")#::Vector{ Vector{FusionRing} }
  product == "TensorProduct" ||
    error("Only tensor product decompositions are defined at the moment.")

  tpd = fr.tensor_product_decompositions
  if tpd !== missing
    [[awc(code) for code in decomp] for decomp in tpd]
  else
    tensor_product_decompositions(fr)
  end
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           tensor_product_decompositions                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
export tensor_product_decompositions

###
"""
    tensor_product_decompositions(r::FusionRing; digits::Int=10)

Return decompositions of `r` as tensor products of fusion rings.

combines:
1. known-ring search;
2. internal subring discovery;
3. Frobenius-Perron dimension signatures.

FPdims are used as  necessary tensor-product invariant:

    FPdim_R simples = products of FPdim_factor simples.

 FPdims alone do not prove a tensor product decomposition.
Every returned decomposition is still verified using `is_equivalent_fusion_ring`.


returns a vector of decompositions, where each decomposition is itself
       a vector of fusion rings:

           [
               [A, B],
               [C, D, E],
               ...
           ]
"""
function tensor_product_decompositions(r::FusionRing; digits::Int = 10)
  decomps = Vector{Vector{FusionRing}}()

  append!(decomps, known_tensor_product_decompositions(r; digits = digits))
  append!(decomps, discover_tensor_product_decompositions(r; digits = digits))

  return unique_decomps!(decomps)
end

#reminder to self (so when I come back to this later I won't forget):
# Mathematical idea:
#   If R is  tensor product of factors A, B, ..., then  simple objects should
#   behave like tuples
#
#       (a, b, ...)
#
#    a is  simple object of A, b is a simple object of B ..etc
#
#   necessary rank condition is:
#
#       rank(R) = rank(A) * rank(B) * ...
#
#    condition is useful but weak. Many unrelated fusion rings can have
#   compatible ranks.
#
#       FPdim_R(a, b) = FPdim_A(a) * FPdim_B(b).
#
#   Tif R equiv A cross B,  multiset of FPdims of R =  multiset of
#    prods FPdim_A(a) * FPdim_B(b)
#
#   () a ranges over simples of A and b ranges over simples of B.)

"""_multiplicative_partitions(n)
    Returns all unordered multiplicative decompositions of n into factors     >= 2.

       ex:
          _multiplicative_partitions(12)
       includes:
           [2, 6], [3, 4], [2, 2, 3], [12]

       use this because if rank(R) = 12, then possible tensor-factor ranks
      could be 2 and 6, or 3 and 4, or 2, 2, and 3.

"""
function multiplicative_partitions(n::Int; minfactor::Int = 2)
  n == 1 && return [Int[]]

  out = Vector{Vector{Int}}()

  function go(rem::Int, start::Int, acc::Vector{Int})
    if rem == 1
      push!(out, copy(acc))
      return nothing
    end

    for d in start:rem
      rem % d == 0 || continue

      push!(acc, d)
      go(rem ÷ d, d, acc)
      pop!(acc)
    end
  end

  go(n, minfactor, Int[])
  return out
end

"""
    cartesian_choices(lists)

Given a vector of candidate lists, returns all ways to choose one element
from each list.

# Example

If

    [[A1, A2], [B1, B2, B3]]

 passed in,  produces

    [A1, B1], [A1, B2], [A1, B3],
    [A2, B1], [A2, B2], [A2, B3]

 use this when  know  possible ranks of the factors and want to
try all known fusion rings with those ranks.
"""
function cartesian_choices(lists::Vector{Vector{FusionRing}})
  isempty(lists) && return Vector{Vector{FusionRing}}()

  out = Vector{Vector{FusionRing}}()
  cur = Vector{FusionRing}(undef, length(lists))

  function go(i::Int)
    if i > length(lists)
      push!(out, copy(cur))
      return nothing
    end

    for x in lists[i]
      cur[i] = x
      go(i + 1)
    end
  end

  go(1)
  return out
end

#I'm pretty sure we already have these lol

function raw_fpdims(fr::FusionRing)
  try
    return fpdims(fr)
  catch
    return frobenius_perron_dimensions(fr)
  end
end

function numeric_fpdim_values(fr::FusionRing)
  try
    return numeric_fpdims(fr)
  catch
    return raw_fpdims(fr)
  end
end

#        Converts an FPdim value to a real Float64.
#
#       This is needed because FPdims may come from OSCAR/QQBar or other exact
#       algebraic number types, not just ordinary Julia floats.
function to_float_real(x)
  try
    return Float64(real(x))
  catch
    return Float64(x)
  end
end

#   fpdim_signature(fr; digits=10)
#       Computes a sorted, rounded multiset of the FPdims of fr.
#
#       Example:
#           [1, sqrt(2), 1]
#
#       becomes approximately:
#           [1.0, 1.0, 1.4142135624]
#
#       Sorting makes the signature independent of the order of simple objects.
#       Rounding avoids tiny numerical differences from eigenvalue computation.
function fpdim_signature(fr::FusionRing; digits::Int = 10)
  ds = _numeric_fpdim_values(fr)
  vals = [_to_float_real(d) for d in ds]
  return sort(round.(vals; digits = digits))
end

#   fpdim_product_signature(factors; digits=10)
#       Computes the FPdim signature that the tensor product of the given factors
#       would have.
#
#       For factors A and B, it computes all products:
#
#           FPdim_A(a) * FPdim_B(b).
#
#       For more than two factors  computes all products across all factors.
#
#        lets us test whether a proposed list of factors could possibly have
#       tensor product equivalent to R.
function fpdim_product_signature(factors::Vector{FusionRing}; digits::Int = 10)
  vals = [1.0]

  for F in factors
    ds = [_to_float_real(d) for d in _numeric_fpdim_values(F)]
    vals = [x * y for x in vals for y in ds]
  end

  return sort(round.(vals; digits = digits))
end

#   fpdim_signatures_match(R, factors; digits=10)
#       Checks the necessary FPdim condition:
#
#           FPdims(R) == product FPdims(factors).
#
#       If  returns false,  factors cannot tensor to R.
#       If it returns true, the factors are only candidates and still need final
#       verification.
function fpdim_signatures_match(
  R::FusionRing, factors::Vector{FusionRing}; digits::Int = 10
)::Bool
  return fpdim_signature(R; digits = digits) ==
         fpdim_product_signature(factors; digits = digits)
end

#   decomp_key(decomp)
#       Builds a key used to remove duplicate decompositions.
#
#       It uses each factor's rank and FPdim signature rather than barcode,
#       because newly constructed or restricted rings may not have reliable
#       AnyonWiki barcodes.
function decomp_key(decomp::Vector{FusionRing})
  return sort([(rank(R), _fpdim_signature(R)) for R in decomp])
end

#   unique_decomps!(decomps)
#       Removes duplicate decompositions in-place using _decomp_key.
#
function unique_decomps!(decomps::Vector{Vector{FusionRing}})
  unique!(decomp_key, decomps)
  return decomps
end

#   candidate_known_factors_by_rank(R)
#       Groups  known fusion rings frl by rank, keeping only rings whose rank
#       could divide rank(R).
#
#       We use this to cheaply restrict the known-ring search. If rank(K) does
#       not divide rank(R), then K cannot be a tensor factor of R.
#
function candidate_known_factors_by_rank(r::FusionRing)
  Rrank = rank(r)

  by_rank = Dict{Int, Vector{FusionRing}}()

  for K in frl
    rk = rank(K)

    rk <= 1 && continue
    rk == Rrank && continue
    Rrank % rk == 0 || continue

    push!(get!(by_rank, rk, FusionRing[]), K)
  end

  return by_rank
end

#   _known_tensor_product_decompositions(R; digits=10)
#       Searches for tensor product decompositions using known rings from frl.
#
#       Steps:
#           1. Factor rank(R) multiplicatively.
#           2. For each possible factor-rank pattern, choose known rings with
#              those ranks.
#           3. Check the FPdim product signature.
#           4. Construct the tensor product using tensor_product(choice).
#           5. Verify actual fusion-rule equivalence with
#              is_equivalent_fusion_ring(R, T).
#
#       This method can find decompositions when the factors already appear in
#       the known-ring database.
function known_tensor_product_decompositions(r::FusionRing; digits::Int = 10)
  Rrank = rank(r)
  Rrank <= 1 && return Vector{Vector{FusionRing}}()

  by_rank = candidate_known_factors_by_rank(r)
  decomps = Vector{Vector{FusionRing}}()

  for parts in multiplicative_partitions(Rrank)
    length(parts) <= 1 && continue
    all(p -> haskey(by_rank, p), parts) || continue

    lists = [by_rank[p] for p in parts]

    for choice in cartesian_choices(lists)

      # Necessary tensor-product invariant:
      # FPdims of R must =  pairwise/product FPdims of the factors.
      fpdim_signatures_match(r, choice; digits = digits) || continue

      T = tensor_product(choice)

      # FPdims are not enough  - actual fusion rules.
      is_equivalent_fusion_ring(r, T) || continue

      push!(decomps, replace_by_known.(choice))
    end
  end

  return unique_decomps!(decomps)
end

#   all_subring_sets_for_factorization(R)
#       Collects candidate subfusion-ring subsets of R.
#
#        includes:
#           [1]                 - trivial unit subring,
#           sub_fusion_ring_subsets(R),
#           collect(1:rank(R))   whole ring.
#
#        decomposition search later ignores trivial and whole-ring factors,
#       but including them here makes the helper complete and reusable.
# Internal subring discovery
function all_subring_sets_for_factorization(fr::FusionRing)
  r = rank(fr)

  sets = Vector{Vector{Int}}()

  push!(sets, [1])
  append!(sets, sub_fusion_ring_subsets(fr))
  push!(sets, collect(1:r))

  unique!(sets)
  sort!(sets; by = S -> (length(S), S))

  return sets
end

#       Tests whether two subring subsets Aset and Bset multiply like independent
#       tensor factors inside R.
#
#       In a clean tensor product R = A cross B, objects from A and B should satisfy:
#
#           (a, 1) ⊗ (1, b) = (a, b),
#
#       which is single simple object.
#
#        for each a in Aset and b in Bset  product a ⊗ b should:
#           1. have exactly one fusion outcome;
#           2. produce a simple object not already produced by another pair;
#           3. collectively cover every simple object of R.
#
#       If these conditions hold,  products form a Cartesian grid of simples.
#       This is strong evidence that Aset and Bset are internal tensor factors.
#
function unique_product_grid(fr::FusionRing, Aset::Vector{Int}, Bset::Vector{Int})
  r = rank(fr)

  grid = Dict{Tuple{Int, Int}, Int}()
  seen = falses(r)

  @inbounds for a in Aset, b in Bset
    outs = fusion_outcomes(fr, a, b)

    # In  clean tensor product, (a,1) ⊗ (1,b) = (a,b),
    # so  product should be exactly one simple object.
    length(outs) == 1 || return nothing

    x = only(outs)

    # Each pair (a,b) should produce  different simple object.
    seen[x] && return nothing

    grid[(a, b)] = x
    seen[x] = true
  end

  # Product Aset ⊗ Bset should cover all simples of fr.
  all(seen) || return nothing

  return grid
end

#   discover_tensor_product_decompositions(R; digits=10)
#       Searches for decompositions using subfusion rings inside R itself.
#
#       Steps:
#           1. Enumerate candidate subring subsets Aset and Bset.
#           2. Require the rank condition:
#
#                  length(Aset) * length(Bset) == rank(R).
#
#           3. Construct restricted fusion rings A and B.
#           4. Check the FPdim product signature.
#           5. Check the Cartesian product grid condition using
#              unique_product_grid.
#           6. Construct tensor_product([A, B]).
#           7. Verify actual fusion-rule equivalence.
#
#        can find decompositions even when the factors are not already
#       registered as known rings, provided they appear internally as subrings of
#       R.
function discover_tensor_product_decompositions(r::FusionRing; digits::Int = 10)
  Rrank = rank(r)
  Rrank <= 1 && return Vector{Vector{FusionRing}}()

  subrings = all_subring_sets_for_factorization(r)
  decomps = Vector{Vector{FusionRing}}()

  for Aset in subrings, Bset in subrings

    # Ignore trivial and whole-ring factors.
    length(Aset) <= 1 && continue
    length(Bset) <= 1 && continue
    length(Aset) == Rrank && continue
    length(Bset) == Rrank && continue

    # Necessary rank condition for R equiv A cross B.
    length(Aset) * length(Bset) == Rrank || continue

    A = restrict_subring(r, copy(Aset); check_closed = true)
    B = restrict_subring(r, copy(Bset); check_closed = true)

    # Necessary FPdim condition for "   "
    fpdim_signatures_match(r, [A, B]; digits = digits) || continue

    #  check whether products a ⊗ b form  unique Cartesian grid.
    grid = _unique_product_grid(r, Aset, Bset)
    grid === nothing && continue

    T = tensor_product([A, B])

    is_equivalent_fusion_ring(r, T) || continue

    push!(decomps, [replace_by_known(A), replace_by_known(B)])
  end

  return unique_decomps!(decomps)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                adjoint_fusion_ring                              ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export adjoint_fusion_ring

function adjoint_fusion_ring(
  ring::FusionRing; represent_by_known = true
)::Tuple{Vector{Int}, FusionRing}
  d(i) = conjugate_element(ring, i)

  el_seen = falses(rank(ring))
  for (i, j, c) in nzsc(ring)
    if j == d(i)
      el_seen[c] = true
    end
  end
  el = findall(el_seen)

  generatedEl = _fusion_closure(ring, el)

  rbk = represent_by_known ? replace_by_known : identity

  if length(generatedEl) == rank(ring)
    return (generatedEl, rbk(ring))
  else
    sr = restrict_subring(ring, generatedEl; check_closed = true)
    return (generatedEl, rbk(sr))
  end
end

# closure of a subset of elements of a fusion ring under fusion
function _fusion_closure(fr::FusionRing, S0::Vector{Int})::Vector{Int}
  r = rank(fr)
  seen = falses(r)
  @inbounds for s in S0
    seen[s] = true
  end
  changed = true
  while changed
    changed = false
    current = findall(seen)
    @inbounds for a in current, b in current
      for c in fusion_outcomes(fr, a, b)
        if !seen[c]
          seen[c] = true
          changed = true
        end
      end
    end
  end
  return findall(seen)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                               upper_central_series                              ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export upper_central_series

function upper_central_series(fr::FusionRing)
  chain = Tuple{Vector{Int}, FusionRing}[]
  push!(chain, (collect(1:rank(fr)), fr))
  while true
    S, adj = adjoint_fusion_ring(last(chain)[2])
    last(chain)[2] === adj && break
    push!(chain, (S, adj))
    length(S) == 1 && break
  end
  return chain
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  is_nilpotent                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export is_nilpotent

function is_nilpotent(r::FusionRing)::Bool
  chain = upper_central_series(r)
  last_set, last_ring = last(chain)
  return length(last_set) == 1
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  adjoint_irreps                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export adjoint_irreps

"""
    adjoint_irreps(ring::FusionRing) -> Vector{Vector{Int}}

Return the partition of the simple objects of `ring` into subsets invariant
under left and right action of the adjoint subring.

- compute  adjoint subring `(subEl, subRing) = adjoint_fusion_ring(ring)`
- for each simple `e`, take  fixed-point closure of `{e}` under:
      X -> Action[(ring, adjoint), X] ∪ Action[X, (ring, adjoint)]
- remove duplicates

More simply, let (S,A) be adjoint subring of R where S is set of simple indices in that subring
For subset X of simples def:
    Left action: S*X = \$⋃_{a∈S, x∈X} Sup(a⊗x)
    Right action: X*S = \$⋃_{x∈X, a∈S} Sup(x⊗a)
    combined : \$\\theta(X) = S*X ∪ X*S\$
So for each simple e I start w/ {e} and apply that until stabilizes:
    X₀ = {e}
    X₁ = θ(X₀)
    X₂ = θ(X₁)
    ...
    Xₙ = θ(Xₙ₋₁) until Xₙ == Xₙ₋₁
    Then obv remove duplicates at the end

Note: the code is slightly different than that in anyonica but ive highlighted/commented the parts that match
"""
function adjoint_irreps(fr::FusionRing)::Vector{Vector{Int}}
  adj_el, adj_ring = adjoint_fusion_ring(fr)

  r = rank(fr)

  # Functions that return all elements obtained by acting with adj_ring on some elements
  # Left action
  function left_action(elements::Vector{Int})::Vector{Int}
    seen = falses(r)
    @inbounds for a in adj_el, b in elements
      for c in fusion_outcomes(fr, a, b)
        seen[c] = true
      end
    end
    return findall(seen)
  end

  # Right action
  function right_action(elements::Vector{Int})::Vector{Int}
    seen = falses(r)
    @inbounds for b in elements, a in adj_el
      for c in fusion_outcomes(fr, b, a)
        seen[c] = true
      end
    end
    return findall(seen)
  end

  # Combined action
  function combined_action(elements::Vector{Int})::Vector{Int}
    return sort!(unique!(vcat(left_action(elements), right_action(elements))))
  end

  # FixedPoint[CombinedAction[pair], [e]]
  function closure(seed::Int)::Vector{Int}
    cur = [seed]
    while true
      nxt = combined_action(cur)
      nxt == cur && return cur
      cur = nxt
    end
  end

  return unique(closure.(collect(1:r)))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                universal_grading                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export universal_grading

function universal_grading(fr::FusionRing)
  irreps = adjoint_irreps(fr)
  n = size(irreps, 1)

  grading = Pair{Int, Int}[]
  @inbounds for a in 1:n
    for x in irreps[a]
      push!(grading, x => a)
    end
  end
  sort!(grading; by = p -> first(p))

  function _cond(l1::Vector{Int}, l2::Vector{Int}, l3::Vector{Int})::Bool
    @inbounds for i in l1, j in l2
      for c in fusion_outcomes(fr, i, j)
        c ∉ l3 && return false
      end
    end
    return true
  end

  mt = zeros(Int, n, n, n)
  @inbounds for a in 1:n, b in 1:n, c in 1:n
    mt[a, b, c] = _cond(irreps[a], irreps[b], irreps[c]) ? 1 : 0
  end

  groupRing = fusion_ring(mt)
  return grading, replace_by_known(groupRing)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   all_gradings                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
export all_gradings

"""
    all_gradings(fr::FusionRing)

Return all gradings of fr

Every grading  obtained from  universal grading by quotienting
universal grading group by a normal subgroup.

    C = ⨁_{g ∈ U} C_g

is  universal grading,  every grading  obtained from U/N,
where N is a normal subgroup of U.
"""

function all_gradings(fr::FusionRing)
  ag = fr.all_gradings
  !ismissing(ag) && return ag

  grading, universal_group = universal_grading(fr)

  universal_blocks = blocks_from_grading_pairs(grading)

  G    = to_group(universal_group)
  mult = cayley_table(universal_group)

  gradings = Vector{Dict{String, Any}}()

  for N in normal_subgroups(G)
    cosets = cosets_as_index_blocks(G, N)

    partition = merge_universal_blocks(universal_blocks, cosets)

    quotient_mult = quotient_mult_table(mult, cosets)
    quotient_group_ring = group_ring_from_cayley_table(quotient_mult)

    Q, proj = quo(G, N)

    push!(
      gradings,
      Dict{String, Any}(
        "partition" => partition,
        "group_ring" => quotient_group_ring,
        "normal_subgroup" => N,
        "quotient_group" => Q,
        "quotient_projection" => proj,
        "cosets" => cosets,
      ),
    )
  end

  sort!(gradings; by = g -> -length(g["partition"]))

  return gradings
end

function blocks_from_grading_pairs(grading::Vector{Pair{Int, Int}})
  ngrades = maximum(last.(grading))
  blocks = [Int[] for _ in 1:ngrades]

  for p in grading
    simple = first(p)
    grade = last(p)
    push!(blocks[grade], simple)
  end

  foreach(sort!, blocks)
  return blocks
end

function merge_universal_blocks(
  universal_blocks::Vector{Vector{Int}}, cosets::Vector{Vector{Int}}
)
  partition = Vector{Vector{Int}}()

  for C in cosets
    block = Int[]

    for g in C
      append!(block, universal_blocks[g])
    end

    push!(partition, sort(unique(block)))
  end

  sort!(partition; by = B -> minimum(B))

  return partition
end

function quotient_mult_table(mult::AbstractMatrix{<:Integer}, cosets::Vector{Vector{Int}})
  q = length(cosets)

  element_to_coset = Dict{Int, Int}()

  for (i, C) in pairs(cosets)
    for g in C
      element_to_coset[g] = i
    end
  end

  quotient_mult = zeros(Int, q, q)

  for A in 1:q
    for B in 1:q
      a = cosets[A][1]
      b = cosets[B][1]
      c = mult[a, b]

      quotient_mult[A, B] = element_to_coset[c]
    end
  end

  return quotient_mult
end

function group_ring_from_cayley_table(mult::AbstractMatrix{<:Integer})
  q = size(mult, 1)

  mt = zeros(Int, q, q, q)

  for a in 1:q
    for b in 1:q
      c = mult[a, b]
      mt[a, b, c] = 1
    end
  end

  return replace_by_known(fusion_ring(mt))
end
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    commutator                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export commutator

""" commutator(fr, sub) -> (els, subring)

Centralizer-style commutator of subring `sub = (S, RS)` inside `fr`:
return all simples `i` with `FusionOutcomes(i ⊗ i*) ⊆ S`
"""

function commutator(fr::FusionRing, sub::Tuple{Vector{Int}, FusionRing})
  is_commutative(fr) || error("commutator: ring must be commutative")

  subEls, _ = sub
  Sset = Set(subEls)

  in_sub(i::Int)::Bool = begin
    di = conjugate_element(fr, i)
    outs = fusion_outcomes(fr, i, di)
    all(in(Sset), outs)
  end

  els0 = [i for i in 1:rank(fr) if in_sub(i)]
  isempty(els0) && (els0 = [1])

  els = _fusion_closure(fr, els0)
  return els, restrict_subring(fr, els; check_closed = true)
end

function commutator(fr::FusionRing)
  return commutator(fr, (collect(1:rank(fr)), fr))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    characters                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export characters

function characters(ring::FusionRing; use_numerics = true)
  if !(ring.characters === missing)
    return from_qqb_id(ring.characters)
  elseif !FusionRings.is_commutative(ring)
    error(
      "Calculation of characters for non-commutative fusion ring is not implemented yet."
    )
  elseif rank(ring) == 1
    return [qqbar(1)]
  else
    mt   = multiplication_table(ring)
    r    = rank(ring)
    mats = [mt[i, :, :] for i in 1:r]
    qqb  = algebraic_closure(QQ)

    sort_mat(mat) = sortslices(mat, dims = 1, by = char_sort_crit)

    if use_numerics
      # find numeric diagonalizing matrix
      V = (normalize_first_col ∘ numeric_diagonalizing_matrix)(mats)
      # find symbolic evals
      evals = union(vcat([eigenvalues(matrix(qqb, m)) for m in mats] ...))
      # identify evals in V, sort mat, and convert to matrix of qqb's
      chars = replace_by_known(evals).(V)
    else
      chars = (normalize_first_col ∘ diagonalizing_matrix)(mats)
    end
    matrix(qqb, sort_mat(chars))
  end
end

function numeric_diagonalizing_matrix(mats, tries::Int = 64, tol::Real = 1e-12)
  r = length(mats)

  function is_diagonalizing_matrix(mat, mats)
    imat = inv(mat)
    for i in 1:r
      D = imat * mats[i] * mat
      @inbounds for j in 1:r
        ;
        D[j, j] = 0.0;
      end #set diagonal el = 0

      if norm(D) > tol # D should be 0 mat now
        return false
      else
        continue
      end
    end
    return true
  end

  for i in 1:tries
    coeffs = rand(Float64, r)
    M = zeros(Float64, r, r)
    @inbounds for k in 1:r
      M .+= coeffs[k] .* mats[k]
    end

    V = Matrix(eigen(M).vectors)    # symmetric not guaranteed; generic eigen

    # need to take inv sicne Oscar's (and our) definition of diagonalizing matrix
    # is inverse of that of LinearAlgebra package
    is_diagonalizing_matrix(V, mats) && return inv(V)

    continue
  end

  return error(
    "No diagonalizing matrix found in " *
    string(tries) *
    " tries. Try setting tries to a higher value or set a different seed for random generator.",
  )
end

function diagonalizing_matrix(mats)
  qqbmats = [matrix(algebraic_closure(QQ), m) for m in mats]

  function is_diagonalizing_matrix(mat, mats)
    invmat = inv(mat)
    for m in mats
      if !is_diagonal(mat * m * invmat)
        return false
      else
        continue
      end
    end
    return true
  end

  proposed_mat = qqbmats[1]

  r = first(size(first(mats)))

  diagq = false;
  upi = 4;
  upj = 4

  while !diagq
    upi += 1
    upj += 1

    # Take random linear rational combination of matrices in mats
    rvec = rand(unique([i//j for i in 1:upi, j in 1:upj]), r)
    sgnvec = rand([-1 1], r)
    combinedmat = sgnvec[1] * rvec[1] * qqbmats[1]
    for i in 1:r
      combinedmat += sgnvec[i] * rvec[i] * qqbmats[i]
    end

    # Find diagonalizing matrix
    proposed_mat = reduce(vcat, (collect ∘ values ∘ eigenspaces)(combinedmat))

    # Check whether it works
    diagq = is_diagonalizing_matrix(proposed_mat, qqbmats)
  end

  return proposed_mat
end

# Sort criterion for characters
function char_sort_crit(v)
  RR = ArbField(64);
  CC = AcbField(64);
  conv(x) = convert(Float64, x)
  # Abs values of elements of v
  absval(vec) = conv.(RR.(abs2.(vec)))
  # Angles of elements of v
  angl(vec) = conv.(real.(log.(CC.(vec)) ./ CC(2 * pi * im)))

  return (- Int(all(isreal.(v))), angl(v), absval(v))
end

function normalize_first_col(mat)
  m, n = size(mat)
  return [mat[i, j] / mat[i, 1] for i in 1:m, j in 1:n]
end

function is_diagonalizing_matrix(mat, ring::FusionRing)
  mt   = FusionRings.multiplication_table(ring)
  r    = FusionRings.rank(ring)
  mats = [matrix(qqb, mt[i, :, :]) for i in 1:r]
  mat  = matrix(qqb, mat)

  return all(is_diagonal(mat * m * inv(mat)) for m in mats)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                formal_codegrees                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export formal_codegrees

function formal_codegrees(fr::FusionRing)::Vector{QQBarFieldElem}
  chars = characters(fr)
  d = conjugate_element(fr)

  return sum(chars[i, :]' * chars[d(i), :] for i in 1:r)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                               numeric_characters                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export numeric_characters
"""
    numeric_characters(R::FusionRing; tries::Int=64, tol=1e-12) -> (C, V)

#Return the **character table** `C::Matrix{ComplexF64}` of a **commutative** fusion ring `R`.

#Algorithm:
#1. Form a random real combination `M = ∑_k c_k N_k`.
#2. Eigen-decompose `M = V Λ V⁻¹`.
#3. Verify that every `V⁻¹ N_i V` is (numerically) diagonal. If not, retry.
#Normalize and sort V

Throws if no common eigenbasis is found after `tries` attempts.
"""
function numeric_characters(ring, tries::Int = 64, tol = 1e-12)
  if !(ring.numeric_characters === missing)
    return ring.numeric_characters
  elseif !FusionRings.is_commutative(ring)
    error(
      "Calculation of characters for non-commutative fusion ring is not implemented yet."
    )
  elseif rank(ring) == 1
    return [qqbar(1)]
  else
    mt   = multiplication_table(ring)
    r    = rank(ring)
    mats = [mt[i, :, :] for i in 1:r]
    V    = (normalize_first_col ∘ numeric_diagonalizing_matrix)(mats)

    sortslices(V; dims = 1, by = num_char_sort_crit)
  end
end

# Numeric sort criterion for characters
function num_char_sort_crit(v)
  # for normalizing demands
  norm(vec) = sqrt(sum(abs2.(vec)))
  # Check whether all elements are real
  are_real(vec) = (Int ∘ all)(abs(imag(x)) < 1e-14 for x in vec ./ norm(vec))
  # Abs values of elements of v
  absval(vec) = abs2.(vec)
  # Angles of elements of v
  angl(vec) = angle.(vec)

  return (- are_real(v), angl(v), absval(v))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              projective_SL2Z_reps                               ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: function not yet ready for export: not even sure I want to keep it...
#export projective_SL_2_ZZ_reps

function projective_SL_2_ZZ_reps(fr::FusionRing)
  md = fr.projective_SL2Z_reps
  if md !== missing
    return md
  else
    error("No data available and calculation not implemented yet")
  end
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  is_commutative                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export is_commutative

function is_commutative(fr::FusionRing)
  N = multiplication_table(fr)
  r = rank(fr)
  for a in 1:r, b in 1:r, c in 1:r
    N[a, b, c] == N[b, a, c] || return false
  end
  return true
end

"""
  #conjugate_element(fr, a) -> Int

  #Return the integer index of the dual (conjugate) simple object of `a`.
  #Accepts an integer index, a `String`, or a `Symbol`.

#Rationale: internal computations (e.g. composing with other index-based
#operations) are simpler when the result is an index rather than a label.
#Use `conjugate_label` if you need the string form.
"""

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            categories_with_properties                           ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export categories_with_properties

function categories_with_properties(fr::FusionRing)
  return fr.has_categories_with_props
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                is_categorifiable                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
export is_categorifiable

function is_categorifiable(fr::FusionRing)
  return fr.categorifiable
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                 which_injection                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export which_injection
"""
    #which_injection(subring, ring) -> Dict{Int,Int} or nothing

#Find an injection of `subring` into `ring` by:
#1) enumerating fusion-closed subsets `S ⊂ ring` of size rank(subring) containing 1,
#2) comparing multiplication tables up to permutation.
"""

function which_injection(subring::FusionRing, ring::FusionRing)
  rs = rank(subring)
  rr = rank(ring)
  rs > rr && return nothing

  Nbig = multiplication_table(ring)
  Nsub = multiplication_table(subring)

  for S in _internal_closed_subsets(ring, rs)
    Nres = @views Nbig[S, S, S]
    perm = _permutation_vector_equiv(Nsub, Nres)
    perm === nothing && continue

    inj = Dict{Int, Int}()
    @inbounds for i in 1:rs
      inj[i] = S[perm[i]]
    end
    return inj
  end
  return nothing
end

#Added: from updates/commutator
# -  compute the invariant:
#       k(i) = |{ c : N[i,i,c] > 0 }|
# - We partition indices by this invariant.
# - Groups are sorted deterministically (increasing k, then index).

"""
    _internal_closed_subsets(fr, k) -> Vector{Vector{Int}}

Return all fusion-closed subsets of size `k` containing the unit `1`,
generated as `S = [1; T]` where `T` ranges over (k-1)-subsets of `2:r`.
"""
function _internal_closed_subsets(fr::FusionRing, k::Int)::Vector{Vector{Int}}
  r = rank(fr)
  (k <= 0 || k > r) && return Vector{Vector{Int}}()

  if k == 1
    return _internal_multiplication(fr, [1]) ? [[1]] : Vector{Vector{Int}}()
  end

  candidates = Vector{Vector{Int}}()
  for T in combinations(collect(2:r), k - 1)
    S = vcat(1, collect(T))
    _internal_multiplication(fr, S) && push!(candidates, S)
  end
  return candidates
end

# checks whether multiplication of elements in S is internal in fusion ring fr
function _internal_multiplication(fr::FusionRing, S::Vector{Int})::Bool
  Sset = Set(S)
  @inbounds for i in S, j in S
    for c in fusion_outcomes(fr, i, j)
      c in Sset || return false
    end
  end
  return true
end

#Added from: updates/automorphisms_which_injections
# generate all k-subsets of 1:n without external deps - could not find corresponding combinatorics function
# so can be replaced by function once found
function _k_subsets(n::Int, k::Int)
  out = Vector{Vector{Int}}()
  buf = Vector{Int}(undef, k)
  function go(start::Int, depth::Int)
    if depth > k
      push!(out, copy(buf));
      return nothing
    end
    # ensure enough remaining
    last = n - (k - depth)
    for v in start:last
      buf[depth] = v
      go(v + 1, depth + 1)
    end
  end
  k == 0 && return [Int[]]
  (k < 0 || k > n) && return out
  go(1, 1)
  return out
end
