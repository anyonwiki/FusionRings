
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                          various helper functions                               ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# TODO: function is not used except in tests...
"Return the fusion matrix (left multiplication by `a`)."
function fusion_matrix(fr::FusionRing, a::Int)::Matrix{Int}
    @views multiplication_table(fr)[a, :, :]
end

"Structure constant N[a,b,c]."
# TODO: function is not used except in tests...
function fusion_coeff(fr::FusionRing, a::Int, b::Int, c::Int)::Int
    multiplication_table(fr)[a,b,c]
end

"""
    fusion_product(fr, a, b) -> Dict{Int,Int}

Return the decomposition of `a ⊗ b` as a multiplicity dictionary
`Dict{simple_index => multiplicity}`.
"""
function fusion_product(fr::FusionRing, a::Int, b::Int)
    N = @views multiplication_table(fr)[a,b,:]
    out = Dict{Int,Int}()
    @inbounds for (c,m) in enumerate(N)
        m==0 && continue
        out[c] = m
    end
    out
end

"Return vector of simple indices with positive multiplicity in `a × b`."
function fusion_outcomes(fr::FusionRing, a::Int, b::Int)::Vector{Int}
    [c for (c,m) in fusion_product(fr,a,b) if m>0]
end

"Ordered list form of `a × b`."
# TODO: function is not used except in tests...
function decompose(fr::FusionRing, a::Int, b::Int) 
    [ (k,v) for (k,v) in fusion_product(fr,a,b) ]
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                replace_by_known                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export replace_by_known

# keep_order=true permutes the elements of the known ring to match the
# elements of the given ring.
# safe_return=true makes the function return the original ring in case no built-in 
# ring is found

function replace_by_known( fr::FusionRing; keep_order=true, safe_return=true )
    r = rank(fr)
    m = multiplicity(fr)
    n = nnsd(fr)

    function possibly_equivalent( R::FusionRing )::Bool 
        fc = anyonwiki_code(R)
        fc[1] == r && fc[2] == m && fc[3] == n && return true
    end

    proposals = filter( possibly_equivalent, frl )

    if isempty(proposals) 
        safe_return && return fr
        return missing
    end

    # find the permutation
    ring_perm_couples = map( ring -> ( ring, which_permutation( ring, fr ) ), proposals )

    target_ind = findfirst( tuple -> tuple[2] !== nothing, ring_perm_couples )

    if target_ind === nothing 
        safe_return && return fr
        return missing
    end

    target = ring_perm_couples[target_ind]

    if keep_order
        permute( target[2][1], target[1] )
    else
        return target
    end

end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                restrict_subring                                 ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export restrict_subring

# restrict ring to subindices S (Anyonica's MT[ring][[el,el,el]])
function restrict_subring(fr::FusionRing, S::Vector{Int}; check_closed::Bool = true)::FusionRing
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
    fusion_ring( Nsub ) 
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                    permute                                      ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: implement permute for permutations from OSCAR
export permute

"""permute(r, perm) – return a new `FusionRing` with all data
    permuted by `perm`.  `perm[1]` **must** equal 1 to keep the vacuum first."""
function permute( perm::Vector{Int}, r::FusionRing)::FusionRing
    n = rank(r)
    n == length(perm)      || throw(ArgumentError("perm length ≠ rank"))
    sort(perm) == collect(1:n) || throw(ArgumentError("perm must be a true permutation"))
    perm[1] == 1 || throw(ArgumentError("vacuum must stay at index 1"))

    # Core table
    mt_new = permute_mult_tab(multiplication_table(r), perm)

    # Metadata that needs re‑ordering (guard against `missing`)
    el_names = labels(r)[perm]
    tex_names = length(r.texnames) == n ? r.texnames[perm] : r.texnames
    fpdims   = r.frobenius_perron_dimensions === missing ?
               missing : (length(r.frobenius_perron_dimensions) == n ? r.frobenius_perron_dimensions[perm] : r.frobenius_perron_dimensions)
    chars    = r.characters === missing ?
               missing : (ndims(r.characters) == 2 && size(r.characters, 2) == n ? r.characters[:, perm] : r.characters)

    return fusion_ring(mt_new;                       # core data
        names         = r.names,
        texnames      = tex_names,
        labels        = el_names,
        barcode       = r.barcode,
        anyonwiki_code = r.anyonwiki_code,
        sub_fusion_rings = r.sub_fusion_rings,
        frobenius_perron_dimensions = fpdims,
        characters    = chars
    )
end

"""
    permute_mult_tab(N, p)

Apply permutation `p` (fixing 1) to all three indices of `N`.
"""
function permute_mult_tab(N::Array{Int,3}, p::Vector{Int})
    p[1]==1 || error("Permutation must fix the unit at index 1")
    r = size(N,1)
    M = fill(0, r, r, r)
    @inbounds for a in 1:r, b in 1:r, c in 1:r
        M[p[a], p[b], p[c]] = N[a,b,c]
    end
    M
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                     sort                                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#TODO: figure out which one is best: string or symbol as named input 
#TODO: make sort function have secondary sort method based on e.g. grouping dual 
# elements together

export sort

function sort( fr::FusionRing; by="fpdims", order::Symbol = :increasing )
    if by == "fpdims"
        return permute( perm_vec_qd(fr,order=order), fr )
    elseif by == "sd_conj"
        return permute( perm_vec_sd_conj(fr,order=order), fr )
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

"""perm_vec_sd_conj(r; order = :increasing) – self‑duals first, then conjugate
    pairs, each block ordered by FP‑dimension."""
function perm_vec_sd_conj(r::FusionRing; order::Symbol = :increasing)::Vector{Int}
    n  = rank(r)
    conj = conjugate_element(r)
    qd = fpdims(r)

    self_dual = [i for i in 2:n if conj(i) == i]
    sort!(self_dual; by = i -> qd[i], rev = (order == :decreasing))

    paired   = Set(self_dual)
    pairs    = Tuple{Int,Int}[]

    for i in 2:n
        i in paired && continue
        j = conj(i)
        i == j && continue

        a, b = i, j
        if (order == :increasing && qd[a] > qd[b]) || (order == :decreasing && qd[a] < qd[b])
            a, b = b, a
        end

        push!(pairs, (a, b))
        push!(paired, a)
        push!(paired, b)
    end

    sort!(pairs; by = p -> qd[p[1]], rev = (order == :decreasing))
    conjlist = reduce(vcat, ([p[1], p[2]] for p in pairs); init = Int[])

    vcat(1, self_dual, conjlist)
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                tensor_product                                   ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

export tensor_product

function tensor_product(r1::FusionRing, r2::FusionRing)::FusionRing
    m, n   = rank(r1), rank(r2)
    mt1, mt2 = multiplication_table(r1), multiplication_table(r2)

    mt = zeros(Int, m*n, m*n, m*n)
    @inbounds for a in 1:m, α in 1:n, b in 1:m, β in 1:n, c in 1:m, γ in 1:n
        i = (a-1)*n + α
        j = (b-1)*n + β
        k = (c-1)*n + γ
        mt[i,j,k] = mt1[a,b,c] * mt2[α,β,γ]
    end

    # Assemble element names
    elnames = [ string(e1, "⊗", e2) for e1 in labels(r1) for e2 in labels(r2) ]

    names_tp = (isempty(names(r1)) || isempty(names(r2))) ? String[] :
        [string(names(r1)[1], "⊗", names(r2)[1])]

    fpdims_new = try
        [d1 * d2 for d1 in fpdims(r1) for d2 in fpdims(r2)]
    catch
        missing
    end

    return fusion_ring(
        mt; 
        names = names_tp, 
        labels = elnames, 
        frobenius_perron_dimensions = fpdims_new
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

function to_group( fr::FusionRing )
  ct = cayley_table(fr)
  gm( a, b ) = ct[a,b]
  permutation_group( generic_group(1:rank(fr), gm )[1] )
end



# TODO: implement bicrossed product. @Szagha02: not a priority 
# @gvercley will do this at some point