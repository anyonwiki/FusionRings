
#function change_fusion_ring_property(r::FusionRing, dict)

#end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            multiplication_table                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export multiplication_table

function multiplication_table(r::FusionRing)::Array{Int,3}
  return r.multiplication_table
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    rank                                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export rank

function rank(r::FusionRing)::Int
  size(multiplication_table(r),1)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    names                                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export names 

function names(r::FusionRing)::Array{String,1}
  return r.names
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  tex_names                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export tex_names

function tex_names(r::FusionRing)::Array{String,1}
  return r.texnames
end

 
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    labels                                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export labels

function labels(r::FusionRing)::Array{String,1}
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

export multiplicity

function multiplicity(r::FusionRing)::Int
  maximum(multiplication_table(r))
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      mult                                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export mult

mult = multiplicity

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           nonzero_structure_constants                           ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export nonzero_structure_constants

function nonzero_structure_constants(r::FusionRing)::Vector{Tuple{Int64, Int64, Int64}}
  mt = multiplication_table(r)
  Tuple.( findall( x -> x > 0, mt ) ) 
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

function frobenius_perron_dimensions(r::FusionRing)::Vector{QQBarFieldElem}
  stored_dims = r.frobenius_perron_dimensions 
  if stored_dims===missing
    mt = multiplication_table(r)
    multmats = [ matrix( ZZ, mt[i,:,:] ) for i in 1:rank(r) ]
    return [ first(eigenvalues(QQBar, A)) for A in multmats ]
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
  return sum( fpdims(r).^2 )
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
    v ./ v[1]
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   numeric_fpdim                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export numeric_fpdim

numeric_fpdim(fr::FusionRing) = sum(x->x*x, numeric_fpdims(fr))

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           num_self_dual_non_self_dual                           ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export num_self_dual_non_self_dual

function num_self_dual_non_self_dual(r::FusionRing)::Array{Int,1}
  sd  = count( x -> x == 1, diag( conjugation_matrix(r) ) )
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
  first( nsdnsd(r) )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                      nsd                                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export nsd

nsd = num_self_dual

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                num_non_self_dual                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export num_non_self_dual 

function num_non_self_dual(r::FusionRing)::Int
  last( nsdnsd(r) )
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
  return sum( multiplication_table(r) ) == rank(r)^2
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   cayley_table                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export cayley_table

function cayley_table( fr::FusionRing )
  !is_group_ring(fr) && message("Ring must be group ring")

  mt = multiplication_table(fr)
  r  = rank(fr)
  return [ findfirst( ==(1), mt[a,b,:] ) for a in 1:r, b in 1:r ]
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
    findfirst(==(1), C[a, :])
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                 conjugate_pairs                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export conjugate_pairs

function conjugate_pairs( fr::FusionRing )::Vector{Vector{Int}}
    d  = conjugate_element(fr);
    us = unique ∘ sort
    unique([ us( [ a, d(a) ] ) for a in 1:rank(fr) ])
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                 anyonwiki_code                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export anyonwiki_code

function anyonwiki_code(r::FusionRing)::Union{Array{Int,1},Missing}
  return r.anyonwiki_code
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    barcode                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export barcode 

function barcode(r::FusionRing)
  return r.barcode
end

function mult_tab_code(mat::Array{Int,2},mult::Int)::Int
    error("mult_tab_code not implemented yet")
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                               sub_fusion_rings                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export sub_fusion_rings

function sub_fusion_rings(r::FusionRing)
    dictvec = r.sub_fusion_rings
    if dictvec !== missing
        [
            Dict(
                "injection"   => dict["injection"],
                "fusion_ring" => awc(dict["anyonwiki_code"])
            )
            for dict in dictvec
        ]
    else
        subsets = sub_fusion_ring_subsets
        [
            Dict( 
                "injection"   => s,
                "fusion_ring" => replace_by_known(fusion_ring(mt[s,s,s]))
            )
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
    np    = size(pairs,1)

    # flatten conjugate_pairs back to 1D vector
    function flatten(l::Vector{Vector{Int}})
        res = Int[]
        for vec ∈ l, el ∈ vec 
            push!(res,el)
        end
        return res
    end
    
    for k ∈ 1:(np-1), subset in combinations(pairs, k)
        S = vcat( [1], flatten(collect(subset)) )
        if is_sub_fusion_ring(fr, S )
            push!(result, S)
        end
    end
    return result
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              is_sub_fusion_ring                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export is_sub_fusion_ring

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

    mt = multiplication_table(fr)[S,S,S]

    check_struct_const(mt) && 
    check_mt_dims(mt) &&
    check_unit(mt) && 
    check_inverse(mt) &&
    check_associativity(mt)
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            is_equivalent_fusion_ring                            ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export is_equivalent_fusion_ring

function is_equivalent_fusion_ring(ring1::FusionRing,ring2::FusionRing)::Bool
   which_permutation(ring1,ring2) !== nothing
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                which_permutation                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export which_permutation

"""
    which_permutation(fr1, fr2, all=false, short_circuit=true)::Union{Nothing,Vector{Int64}} 

    Returns a list with a permutation such that `permute( perm, fr1 ) == fr2` or nothing if none exist.
    If all=true, returns all permutations such that `permute( perm, fr1 ) == fr2` or nothing if none exist.
    If short_circuit=false, doesn't compute invariants to determine whether rings are  
    not equivalent. This saves time if it is known a priori that the rings have to be equivalent.
"""
function which_permutation(fr1::FusionRing, fr2::FusionRing; 
    all=false, 
    short_circuit = true)::Union{Nothing,Vector{Vector{Int64}}}

    r = rank(fr1)

    # if ranks are different we can't compare so necessary to short_circuit
    # since computation of number selfdual & non-selfdual elements is as 
    # fast and gives stronger invariant, we use that one instead 
    nsdnsd(fr1) != nsdnsd(fr2) && return nothing 


    # check whether rings have same multiplicity and number of selfdual elements
    if short_circuit && ( mult(fr1) != mult(fr2) || fpdim(fr1) != fpdim(fr2) )
        return nothing
    end

    # check whether number of fusion outcomes per multiplicity 
    # is the same for each element of fr1 and fr2 
    m1 = multiplication_table(fr1)
    m2 = multiplication_table(fr2)

    grp1 = diag_channel_count(m1)
    grp2 = diag_channel_count(m2)

    short_circuit && sort(grp1) ≠ sort(grp2) && return nothing

    # rings must have same fpdims
    fpd1 = fpdims(fr1)
    fpd2 = fpdims(fr2)

    short_circuit && sort(fpd1) ≠ sort(fpd2) && return nothing

    # construct invariants for all elements of both rings
    # these are couples of their fusion outcome count and 
    # their fpdims

    inv1 = collect( zip( grp1, fpd1 ) )
    inv2 = collect( zip( grp2, fpd2 ) )

    #sort invA and invB such that equal elements are next to each other.
    σ1 = sortperm(inv1)
    σ2 = sortperm(inv2)

    sm1 = m1[ σ1, σ1, σ1 ]
    sm2 = m2[ σ2, σ2, σ2 ]

    sorted_inv1 = inv1[σ1] 

    S, _ = symmetries(sorted_inv1, sorted = true)
    
    if !all  # only need first permutation
        for σ ∈ S 
            p = Vector(σ,r)
            if sm1[p,p,p] == sm2
                iσ2 = invperm(σ2)
                return [ iσ2[p[σ1]] ]
            end
        end
    else # want all permutations 
        allperms = Vector{Int}[]

        for σ ∈ S 
            p = Vector(σ,r)
            if sm1[p,p,p] == sm2
                iσ2 = invperm(σ2)
                push!(allperms, iσ2[p[σ1]])
            end
        end

        return allperms
    end
end

# returns the symmetry group S of the sorted vector v together with 
# the permutation σv that sorts v. The symmetries are thus of the form
# inverseperm(σv) ∘ g ∘ σv

function symmetries( v::AbstractVector; sorted = false )::Tuple{PermGroup,Vector{Int64}}
    is_empty(v) && return nothing

    n = size(v,1)

    n == 1 && return symmetric_group(1)

    if sorted 
        return ( _sorted_symmetries( v ), collect(1:n) )
    else
        σv = sortperm(v)
        vs = v[σv]
        return ( _sorted_symmetries(vs), σv )
    end
end

function _sorted_symmetries( v::AbstractVector)::PermGroup
    # compute the sizes n of the individual S_n
    n = size(v,1)

    degrees = Int64[1]

    g_ind = 1
    for i in 2:n
        if v[i] == v[i-1] 
            degrees[g_ind] = degrees[g_ind] + 1
        else
            push!(degrees,1)
            g_ind = g_ind + 1
        end
    end
    
    return inner_direct_product( symmetric_group.( degrees ) )

end


"""
    _diag_channel_groups(N) -> Vector{Vector{Int}}

Partition indices by invariant k(i)=|{c : N[i,i,c]>0}|.

Return groups in deterministic order:
- increasing k
- increasing indices within each group
"""
function diag_channel_count(N::Array{Int,3})::Vector{Tuple{Vector{Int},Vector{Int}}}
    [ tally( N[i,i,:], sort=true ) for i in 1:size(N,1) ]
end

# Apply permutation P on all three indices: A'[i,j,k] = A[P[i],P[j],P[k]]
function _permute_multtab(A::Array{Int,3}, P::Vector{Int})::Array{Int,3}
    r = size(A, 1)
    B = similar(A)
    @inbounds for i in 1:r, j in 1:r, k in 1:r
        B[i,j,k] = A[P[i], P[j], P[k]]
    end
    B
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                            fusion_ring_automorphisms                            ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export fusion_ring_automorphisms
"""
    fusion_ring_automorphisms(fr) -> Vector{Vector{Int}}

Return all permutations `p` whose action on the indices leave the structure constants invariant.
"""

# TODO: test this
function fusion_ring_automorphisms(fr::FusionRing)
   which_permutation( fr, fr, all= true )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                  decompositions                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export decompositions

function decompositions( fr::FusionRing, product="TensorProduct" )#::Vector{ Vector{FusionRing} }
    product == "TensorProduct" ||  error("Only tensor product decompositions are defined at the moment.")

    tpd = fr.tensor_product_decompositions
    if tpd !== missing
        [ [ awc( code ) for code in decomp ] for decomp in tpd ]
    else
        tensor_product_decompositions(fr)
    end
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                           tensor_product_decompositions                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: implement
export tensor_product_decompositions


function _multiplicative_partitions(n::Int; minfactor::Int=2)
    n == 1 && return [Int[]]

    out = Vector{Vector{Int}}()

    function go(rem::Int, start::Int, acc::Vector{Int})
        if rem == 1
            push!(out, copy(acc))
            return
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

function _cartesian_choices(lists::Vector{Vector{FusionRing}})
    isempty(lists) && return Vector{Vector{FusionRing}}()

    out = Vector{Vector{FusionRing}}()
    cur = Vector{FusionRing}(undef, length(lists))

    function go(i::Int)
        if i > length(lists)
            push!(out, copy(cur))
            return
        end

        for x in lists[i]
            cur[i] = x
            go(i + 1)
        end
    end

    go(1)
    return out
end

function _decomp_key(decomp::Vector{FusionRing})
    sort([barcode(R) for R in decomp])
end



# Strategy:
#   1. Known-ring search:
#       Try tensor products of rings from known_rings().
#
#   2. Discovery search:
#        Look for subfusion rings A,B inside R such that:
#            rank(A) * rank(B) = rank(R)
#        and such that every product a ⊗ b, with a ∈ A and b ∈ B,
#        is a single simple object, and these products cover all simples
#        of R exactly once.
#
#   This detects clean tensor decompositions where the two factors appear
#   inside R as A⊗1 and 1⊗B.



function _cartesian_choices(lists::Vector{Vector{FusionRing}})
    isempty(lists) && return Vector{Vector{FusionRing}}()

    out = Vector{Vector{FusionRing}}()
    cur = Vector{FusionRing}(undef, length(lists))

    function go(i::Int)
        if i > length(lists)
            push!(out, copy(cur))
            return
        end

        for x in lists[i]
            cur[i] = x
            go(i + 1)
        end
    end

    go(1)
    return out
end

function _known_tensor_product_decompositions(r::FusionRing)
    Rrank = rank(r)
    Rrank <= 1 && return Vector{Vector{FusionRing}}()

    by_rank = Dict{Int,Vector{FusionRing}}()
    for K in frl
        rk = rank(K)
        rk <= 1 && continue
        rk == Rrank && continue
        Rrank % rk == 0 || continue
        push!(get!(by_rank, rk, FusionRing[]), K)
    end

    decomps = Vector{Vector{FusionRing}}()

    for parts in _multiplicative_partitions(Rrank)
        length(parts) <= 1 && continue
        all(p -> haskey(by_rank, p), parts) || continue

        lists = [by_rank[p] for p in parts]

        for choice in _cartesian_choices(lists)
            T = _tensor_product_many(choice)
            if is_equivalent_fusion_ring(r, T)
                push!(decomps, replace_by_known.(choice))
            end
        end
    end

    unique!( ds -> sort([rank(x) for x in ds]), decomps )
    return decomps
end

function _all_subring_sets_for_factorization(fr::FusionRing)
    r = rank(fr)

    sets = Vector{Vector{Int}}()
    push!(sets, [1])
    append!(sets, sub_fusion_ring_subsets(fr))
    push!(sets, collect(1:r))

    unique!(sets)
    sort!(sets, by = S -> (length(S), S))
    return sets
end

function _unique_product_grid(fr::FusionRing, Aset::Vector{Int}, Bset::Vector{Int})
    r = rank(fr)

    grid = Dict{Tuple{Int,Int},Int}()
    seen = falses(r)

    @inbounds for a in Aset, b in Bset
        outs = fusion_outcomes(fr, a, b)

        # In a clean tensor product, (a,1) ⊗ (1,b) = (a,b),
        # so the product should be exactly one simple object.
        length(outs) == 1 || return nothing

        x = only(outs)

        # Each pair (a,b) should give a different simple object.
        seen[x] && return nothing

        grid[(a,b)] = x
        seen[x] = true
    end

    # The products Aset ⊗ Bset should cover every simple object of fr.
    all(seen) || return nothing

    return grid
end

function _discover_tensor_product_decompositions(r::FusionRing)
    Rrank = rank(r)
    Rrank <= 1 && return Vector{Vector{FusionRing}}()

    subrings = _all_subring_sets_for_factorization(r)
    decomps = Vector{Vector{FusionRing}}()

    for Aset in subrings, Bset in subrings
        # Ignore trivial and whole-ring factors.
        length(Aset) <= 1 && continue
        length(Bset) <= 1 && continue
        length(Aset) == Rrank && continue
        length(Bset) == Rrank && continue

        # Rank condition for a tensor product.
        length(Aset) * length(Bset) == Rrank || continue

        # Check whether products a⊗b form a unique Cartesian grid.
        grid = _unique_product_grid(r, Aset, Bset)
        grid === nothing && continue

        A = restrict_subring(r, copy(Aset); check_closed = true)
        B = restrict_subring(r, copy(Bset); check_closed = true)

        T = tensor_product(A, B)

        # Final verification.
        is_equivalent_fusion_ring(r, T) || continue

        push!(decomps, [replace_by_known(A), replace_by_known(B)])
    end

    unique!(ds -> sort([rank(x) for x in ds]), decomps )
    return decomps
end

"""
    tensor_product_decompositions(r::FusionRing) -> Vector{Vector{FusionRing}}

Return decompositions of `r` as tensor products of fusion rings.

 combines two methods:

1. Known-ring recognition:
   tries tensor products of rings already available in `known_rings()`.

2. Internal factor discovery:
   searches for subfusion rings A and B inside `r` such that
   every product `a ⊗ b`, with `a ∈ A` and `b ∈ B`, is a unique simple object
   and the products cover all simples of `r`.

The second method can discover actual tensor factors even when they are not
already registered as known rings, provided the factors appear as subrings
inside `r`.
"""
function tensor_product_decompositions(r::FusionRing)
    known = _known_tensor_product_decompositions(r)
    discovered = _discover_tensor_product_decompositions(r)

    decomps = vcat(known, discovered)

    unique!(ds -> sort([rank(x) for x in ds]), decomps)
    return decomps
end

 
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                adjoint_fusion_ring                              ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export adjoint_fusion_ring

#Added: from updates/commutator
export adjoint_fusion_ring
function adjoint_fusion_ring(ring::FusionRing)::Tuple{Vector{Int},FusionRing}
    d(i) = conjugate_element(ring, i)

    el_seen = falses(rank(ring))
    for (i, j, c) in nzsc(ring)
        if j == d(i)
            el_seen[c] = true
        end
    end
    el = findall(el_seen)

    generatedEl = _fusion_closure(ring, el)

    if length(generatedEl) == rank(ring)
        return (generatedEl, ring)
    else
        return (generatedEl, restrict_subring(ring, generatedEl; check_closed = true))
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
    findall(seen)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                               upper_central_series                              ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export upper_central_series

#Added: from branch feat/adjoint_upper
function upper_central_series(fr::FusionRing)
    chain = Tuple{Vector{Int},FusionRing}[]
    push!(chain, (collect(1:rank(fr)), fr))
    while true
        S, adj = adjoint_fusion_ring(last(chain)[2])
        last(chain)[2] === adj && break
        push!(chain, (S, adj))
        length(S) == 1 && break
    end
    chain
end

export is_nilpotent

function is_nilpotent(r::FusionRing)::Bool
    chain = upper_central_series(r)
    # Nilpotent ⇔ iterated adjoint subring reaches the trivial (rank-1) subring.
    last_set, last_ring = last(chain)
    return length(last_set) == 1 && rank(last_ring) == 1
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

More simply (this is mostly for me to remember lol): 
let (S,A) be adjoint subring of R where S is set of simple indices in that subring
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
function adjoint_irreps(ring::FusionRing)::Vector{Vector{Int}}
    adRing = adjoint_fusion_ring(ring)   # Tuple{Vector{Int}, FusionRing}

    # Action[{ring, adRing}, elements]
    function _left_action(
        pair::Tuple{FusionRing,Tuple{Vector{Int},FusionRing}},
        elements::Vector{Int}
    )::Vector{Int}
        fr, (subEl, _) = pair
        seen = falses(rank(fr))
        @inbounds for a in subEl, el in elements
            for c in fusion_outcomes(fr, a, el)
                seen[c] = true
            end
        end
        return findall(seen)
    end

    # Action[elements, {ring, adRing}]
    function _right_action(
        elements::Vector{Int},
        pair::Tuple{FusionRing,Tuple{Vector{Int},FusionRing}}
    )::Vector{Int}
        fr, (subEl, _) = pair
        seen = falses(rank(fr))
        @inbounds for el in elements, a in subEl
            for c in fusion_outcomes(fr, el, a)
                seen[c] = true
            end
        end
        return findall(seen)
    end

    # CombinedAction[{ring, adRing}][elements]
    function _combined_action(
        pair::Tuple{FusionRing,Tuple{Vector{Int},FusionRing}},
        elements::Vector{Int}
    )::Vector{Int}
        sort!(unique!(vcat(
            _left_action(pair, elements),
            _right_action(elements, pair)
        )))
    end

    pair = (ring, adRing)

    # FixedPoint[CombinedAction[pair], [e]]
    function _closure(seed::Int)::Vector{Int}
        cur = [seed]
        while true
            nxt = _combined_action(pair, cur)
            nxt == cur && return cur
            cur = nxt
        end
    end

    blocks = Vector{Vector{Int}}()
    for e in 1:rank(ring)
        blk = _closure(e)
        blk in blocks || push!(blocks, blk)
    end

    return blocks
end


#Added: from updates/commutator
#=
Compute `irreps = adjoint_irreps(fr)` (partition of simples).

Create group object with `n = length(irreps)` elements.
`grading` maps each simple `x` to block index `a`.
Multiplication table on the grading group is:
    mt[a,b,c] = 1  iff  FusionOutcomes(i ⊗ j) ⊆ irreps[c]
for all i ∈ irreps[a], j ∈ irreps[b].
=#

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                universal_grading                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export universal_grading

function universal_grading(fr::FusionRing)
    irreps = adjoint_irreps(fr)
    n = length(irreps)

    grading = Pair{Int,Int}[]
    @inbounds for a in 1:n
        for x in irreps[a]
            push!(grading, x => a)
        end
    end
    sort!(grading, by = p -> first(p))

    function _cond(l1::Vector{Int}, l2::Vector{Int}, l3::Vector{Int})::Bool
        S = Set(l3)
        @inbounds for i in l1, j in l2
            for c in fusion_outcomes(fr, i, j)
                c in S || return false
            end
        end
        true
    end

    mt = zeros(Int, n, n, n)
    @inbounds for a in 1:n, b in 1:n, c in 1:n
        mt[a,b,c] = _cond(irreps[a], irreps[b], irreps[c]) ? 1 : 0
    end

    groupRing = fusion_ring(mt; labels = string.(1:n))
    return grading, replace_by_known(groupRing)
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                        UG                                       ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export UG

UG(fr::FusionRing) = universal_grading(fr)

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                   all_gradings                                  ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# TODO all gradings can be obtained as normal subgroups of the universal grading. No 
# need to generate all partitions here. 

export all_gradings

"""
    all_gradings(fr::FusionRing)

Return all gradings of the fusion ring.

Each grading is returned as:
    (partition, group_ring)

where:
- partition = Vector of blocks (Vector{Vector{Int}})
- group_ring = the grading group as a fusion ring

note to self:
grading is: simples g∈G⨆​Cg​
s.t i ∈ C_{a} , j ∈ C_{b} ⇒ i ⊗ j ⊂ C_{a ⋅ b}
"""
function all_gradings(fr::FusionRing)
    irreps = adjoint_irreps(fr)   # base partition
    n = length(irreps)

    # generate all partitions of {1,...,n}
    parts = _set_partitions(n)

    gradings = []

    for P in parts
        # merge irreps according to partition P
        blocks = [
            sort!(reduce(vcat, irreps[I]))
            for I in P
        ]

        # test if this defines a grading
        mt = zeros(Int, length(blocks), length(blocks), length(blocks))

        function valid_product(a, b)
            S = Set(blocks[b])
            for i in blocks[a], j in blocks[b]
                for c in fusion_outcomes(fr, i, j)
                    c in S || return false
                end
            end
            return true
        end

        valid = true

        for a in eachindex(blocks), b in eachindex(blocks)
            found = false
            for c in eachindex(blocks)
                if valid_product(a, c)
                    mt[a,b,c] = 1
                    found = true
                    break
                end
            end
            if !found
                valid = false
                break
            end
        end

        valid || continue

        push!(gradings, (
            blocks,
            fusion_ring(mt; labels = string.(1:length(blocks)))
        ))
    end

    return gradings
end

#=
note due to bell numbers (number of ways to partition set of n els into non-empty disjoint subsets)
this grows exponentially 
for rank 8 - this is 4140 partitions, for rank 10 - 115975 partitions
but i think i can optimize the partition generation by only generating partitions of the set of irreps (which is usually much smaller than rank) 
and then merging the corresponding blocks of simples, 
which should be much more efficient in practice for the fusion rings we care about
=#
function _set_partitions(n::Int)
    if n == 0
        return [[]]
    end

    result = []

    function backtrack(i, current)
        if i > n
            push!(result, deepcopy(current))
            return
        end

        for j in eachindex(current)
            push!(current[j], i)
            backtrack(i + 1, current)
            pop!(current[j])
        end

        push!(current, [i])
        backtrack(i + 1, current)
        pop!(current)
    end

    backtrack(1, Vector{Vector{Int}}())
    return result
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    commutator                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export commutator

""" commutator(fr, sub) -> (els, subring)

Centralizer-style commutator of subring `sub = (S, RS)` inside `fr`:
return all simples `i` with `FusionOutcomes(i ⊗ i*) ⊆ S`
"""

function commutator(fr::FusionRing, sub::Tuple{Vector{Int},FusionRing})
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
    return els, restrict_subring(fr, els; check_closed=true)
end

function commutator(fr::FusionRing)
    return commutator(fr, ( collect( 1:rank(fr) ), fr ) )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    characters                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export characters 

function characters(ring::FusionRing; use_numerics = true )
  if !(ring.characters === missing)
    return from_qqb_id( ring.characters )
  elseif !FusionRings.is_commutative(ring) 
    error("Calculation of characters for non-commutative fusion ring is not implemented yet.")
  elseif rank(ring) == 1
    return [ qqbar(1) ]
  else
    mt   = multiplication_table( ring )
    r    = rank( ring )
    mats = [ mt[i,:,:] for i in 1:r ]
    qqb  = algebraic_closure(QQ)
    
    sort_mat( mat ) = sortslices( mat, dims = 1, by = char_sort_crit )

    if use_numerics
        # find numeric diagonalizing matrix
        V     = ( normalize_first_col ∘ numeric_diagonalizing_matrix )( mats )
        # find symbolic evals
        evals = union( vcat( [ eigenvalues( matrix( qqb, m ) ) for m in mats ] ... ) )
        # identify evals in V, sort mat, and convert to matrix of qqb's
        chars = replace_by_known(evals).(V)
    else
        chars = ( normalize_first_col ∘ diagonalizing_matrix )( mats )
    end
    matrix( qqb, sort_mat( chars ) )
  end
end

function numeric_diagonalizing_matrix( mats, tries::Int=64, tol::Real=1e-12 )
    r = length( mats )

    function is_diagonalizing_matrix( mat, mats )
        imat = inv(mat)
        for i in 1:r
            D = imat * mats[i] * mat
            @inbounds for j in 1:r; D[j,j] = 0.0; end #set diagonal el = 0

            if norm(D) > tol # D should be 0 mat now
                return false
            else
                continue
            end
        end
        return true
    end

    for i in 1:tries
        coeffs = rand( Float64, r )
        M = zeros( Float64, r, r )
        @inbounds for k in 1:r
            M .+= coeffs[k] .* mats[k]
        end


        V  = Matrix(eigen(M).vectors)    # symmetric not guaranteed; generic eigen

        # need to take inv sicne Oscar's (and our) definition of diagonalizing matrix
        # is inverse of that of LinearAlgebra package
        is_diagonalizing_matrix( V, mats ) && return inv( V )

        continue
    end

    error("No diagonalizing matrix found in " * string(tries) *" tries. Try setting tries to a higher value or set a different seed for random generator.")
end

function diagonalizing_matrix( mats )
  qqbmats = [ matrix( algebraic_closure(QQ), m ) for m in mats ]

  function is_diagonalizing_matrix( mat, mats )
    invmat = inv(mat)
    for m in mats
        if !is_diagonal( mat * m * invmat )
            return false
        else
            continue
        end
    end
    return true
  end

  proposed_mat = qqbmats[1]

  r = first( size( first( mats ) ) )
    
  diagq = false; upi = 4; upj = 4

  while !diagq
    upi += 1
    upj += 1

    # Take random linear rational combination of matrices in mats
    rvec    = rand( unique( [ i//j for i ∈ 1:upi, j ∈ 1:upj ] ), r )
    sgnvec  = rand( [ -1 1 ], r )
    combinedmat = sgnvec[1] * rvec[1] * qqbmats[1]
    for i ∈ 1:r
      combinedmat += sgnvec[i] * rvec[i] * qqbmats[i]
    end

    # Find diagonalizing matrix
    proposed_mat = 
      reduce( 
        vcat,
        (collect ∘ values ∘ eigenspaces)( combinedmat )
      )

    # Check whether it works
    diagq = is_diagonalizing_matrix( proposed_mat, qqbmats )
  end

  return proposed_mat 
end

# Sort criterion for characters
function char_sort_crit( v )
	RR = ArbField(64);
	CC = AcbField(64);
	conv(x) = convert(Float64,x)
	# Abs values of elements of v
	absval(vec) = conv.( RR.( abs2.( vec ) ) )
	# Angles of elements of v
	angl(vec) = conv.( real.( log.( CC.( vec ) ) ./ CC( 2 * pi * im ) ) )
			
	( - Int( all( isreal.(v) ) ) , angl(v), absval(v) )
end

function normalize_first_col( mat )
    m, n = size( mat )
    [ mat[i,j] / mat[i,1] for i in 1:m, j in 1:n ]
end

function is_diagonalizing_matrix( mat, ring::FusionRing )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( qqb, mt[ i, :, : ] ) for i ∈ 1:r ]
  mat  = matrix( qqb, mat )

  all( is_diagonal( mat * m * inv(mat) ) for m in mats )
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
function numeric_characters( ring, tries::Int = 64, tol=1e-12 )
  if !(ring.numeric_characters === missing)
    return ring.numeric_characters
  elseif !FusionRings.is_commutative(ring)
    error("Calculation of characters for non-commutative fusion ring is not implemented yet.")
  elseif rank(ring) == 1
    return [ qqbar(1) ]
  else
    mt   = multiplication_table( ring )
    r    = rank( ring )
    mats = [ mt[i,:,:] for i in 1:r ]
    V    = ( normalize_first_col ∘ numeric_diagonalizing_matrix )( mats )

    sortslices( V, dims = 1, by = num_char_sort_crit )
  end
end

# Numeric sort criterion for characters
function num_char_sort_crit( v )
    # for normalizing demands
    norm(vec) = sqrt( sum( abs2.(vec) ) )
    # Check whether all elements are real
    are_real(vec) = ( Int ∘ all )( abs( imag( x ) ) < 1e-14 for x in vec ./ norm(vec) )
	# Abs values of elements of v
	absval(vec) = abs2.( vec )
	# Angles of elements of v
	angl(vec) = angle.( vec )

	( - are_real(v) , angl(v), absval(v) )
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              projective_SL2Z_reps                               ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: function not yet ready for export: not even sure I want to keep it...
#export projective_SL_2_ZZ_reps

function projective_SL_2_ZZ_reps( fr::FusionRing )
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
        N[a,b,c] == N[b,a,c] || return false
    end
    true
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

function categories_with_properties( fr::FusionRing )
    return fr.has_categories_with_props
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                is_categorifiable                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
export is_categorifiable

function is_categorifiable( fr::FusionRing )
    return fr.categorifiable
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                restrict_subring                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export restrict_subring
#Added: from pushed files branch
# restrict ring to subindices S (Anyonica's MT[ring][[el,el,el]])
function restrict_subring(fr::FusionRing, S::Vector{Int}; check_closed::Bool = true)::FusionRing
    sort!(S)
    N = multiplication_table(fr)
    @views Nsub = N[S, S, S]

    if check_closed
        # sanity: S must be fusion-closed
        rmap = zeros(Int, rank(fr))
        @inbounds for (k, v) in enumerate(S)
            rmap[v] = k
        end
        @inbounds for a in S, b in S
            for c in findall(>(0), N[a, b, :])
                rmap[c] != 0 || error("subset not fusion-closed")
            end
        end
    end

    # TODO: might want to conserve as much information as possible, but 
    # it's better to wait until all fields of the FusionRing struct are finalized
    fusion_ring( Nsub )
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

        inj = Dict{Int,Int}()
        @inbounds for i in 1:rs
            inj[i] = S[perm[i]]
        end
        return inj
    end
    nothing
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
    true
end

#Added from: updates/automorphisms_which_injections
# generate all k-subsets of 1:n without external deps - could not find corresponding combinatorics function
# so can be replaced by function once found
function _k_subsets(n::Int, k::Int)
    out = Vector{Vector{Int}}()
    buf = Vector{Int}(undef, k)
    function go(start::Int, depth::Int)
        if depth > k
            push!(out, copy(buf)); return
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
    go(1,1)
    out
end

