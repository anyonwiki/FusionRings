export psu2k_fusion_ring,
  su2k_fusion_ring,
  son2_fusion_ring,
  metaplectic_fusion_ring,
  group_fusion_ring,
  zn_fusion_ring,
  group_rep_fusion_ring,
  HI_fusion_ring,
  TY_fusion_ring,
  son2_fusion_ring

###############################################################################################
# Fusion Rings related to finite Groups
###############################################################################################

# TODO: we can add more info on group rings
export group_fusion_ring

"""group_fusion_ring(g::Group)::FusionRing takes a Group object g (defined in Oscar) and returns the group ring corresponding to g.
group_fusion_ring(ct::Matrix{Int{64}};name=G)::FusionRing returns the group fusion ring with multiplication corresponding to the Cayley table ct.
"""

function group_fusion_ring(grp::Group)::FusionRing
  ct = cayley_table(grp)
  nm = describe(grp)

  return group_fusion_ring(ct; names = [nm], checktable = false)
end

function group_fusion_ring(
  ct::Matrix{Int64}; names::Vector{String} = ["G"], checktable::Bool = true
)::FusionRing
  if checktable
    _is_group_table(ct) ||
      throw(ArgumentError("table must be a group multiplication table."))
  end

  r = first(size(ct))
  ar(i) = _e(i, r)

  mt = [ar(ct[i, j])[k] for i in 1:r, j in 1:r, k in 1:r]

  return fusion_ring(
    mt;
    labels = [bold_integer(i) for i in 1:r],
    names = names,
    frobenius_perron_dimension = ZZ(r),
    frobenius_perron_dimensions = fill(ZZ(1), r),
  )
end

"""zn_fusion_ring(n::Int)::FusionRing returns the group ring corresponding to the cyclic group with `n` elements.
"""
# TODO: add missing information
function zn_fusion_ring(n::Int)::FusionRing
  mt = fill(0, n, n, n)

  for i in 0:(n - 1), j in 0:(n - 1)
    k = mod(i + j, n)
    mt[i + 1, j + 1, k + 1] = 1
  end

  return fusion_ring(mt; names = ["ℤ" * subscript_integer(n)], labels = string.(0:(n - 1)))
end

function group_rep_fusion_ring(g)
  m         = character_matrix(g)
  invm      = inv(m)
  r         = size(m, 1)
  toint64   = Int64 ∘ ZZ
  χ(a::Int) = m[a, :]

  mt = [toint64(((χ(a) .* χ(b)) * invm)[c]) for a in 1:r, b in 1:r, c in 1:r]

  nms = describe(g)

  return fusion_ring(
    mt; names = ["Rep("*nms*")"], labels = ["χ"*subscript_integer(i) for i in 1:r]
  )
end

# returns the matrix of character values for group g over 
# the field of cyclotomics
function character_matrix(g)
  ch = character_table(g)
  r  = size(collect(ch), 1)

  # convert character table to matrix over cyclotomics
  m = matrix(QQab, [ch[i, j] for i in 1:r, j in 1:r])

  # if trivial char is not first element, set it as first
  triv_row = fill(QQab(1), r)
  ind_triv = findfirst(i -> m[i, :] == triv_row, 1:r)

  if ind_triv != 1
    m[ind_triv, :], m[1, :] = m[1, :], m[ind_triv, :]
  end

  return m
end

#Haagerup–Izumi (HI) and Tambara–Yamagami (TY) fusion rings
"""
    HI_fusion_ring(tab; names=String[]) -> FusionRing

Build the Haagerup–Izumi fusion ring from a *symmetric* group multiplication table `tab`.

Rank is 2n. Objects are:
- 1..n   : group elements
- n+1..2n: "rho*g" sector (s X_g), indexed by g=1..n as n+g.
"""
function HI_fusion_ring(g::Group)::FusionRing
  ct    = cayley_table(g);
  names = ["HI(" * describe(g) * ")"]

  return HI_fusion_ring(ct; names = names, checktable = false)
end

function HI_fusion_ring(
  tab::Matrix{Int64}; names::Vector{String} = ["HI(G)"], checktable::Bool = true
)
  if checktable
    _is_group_table(tab) ||
      throw(ArgumentError("FusionRingHI: table must be a group multiplication table."))
  end
  tab == tab' ||
    throw(ArgumentError("FusionRingHI: multiplication table must be symmetric."))

  n = size(tab, 1)
  r = 2*n
  inv = _inverse_vector(tab)

  mats = Matrix{Int}[]

  # For i in 1..2n build N_i as in  Mathematica Which cases.
  @inbounds for i in 1:r
    Ni = zeros(Int, r, r)
    for j in 1:r
      if i <= n && j <= n
        k = tab[i, j]
        Ni[j, k] += 1
      elseif i <= n && j > n
        # k == n + tab[[i, j-n]]
        k = n + tab[i, j - n]
        Ni[j, k] += 1

      elseif i > n && j <= n
        # k == n + tab[[ inv[[j]], i-n ]]
        k = n + tab[inv[j], i - n]
        Ni[j, k] += 1

      else
        # i>n && j>n:
        # If[ k == tab[[ i-n, inv[[j-n]] ]] || k > n, 1, 0 ]
        # => all "rho-sector" (k>n) appear with multiplicity 1,
        #    plus exactly one group element tab[i-n, inv[j-n]].
        k0 = tab[i - n, inv[j - n]]
        Ni[j, k0] += 1
        for k in (n + 1):r
          Ni[j, k] += 1
        end
      end
    end
    push!(mats, Ni)
  end

  mt = _mats_to_mt(mats)

  # Labels: group elements then rho-sector
  labels = [string(i) for i in 1:n]
  append!(labels, ["ρ"*subscript_integer(i) for i in 1:n])

  return fusion_ring(mt; names = names, labels = labels)
end

# Helper function for HI_fusion_ring
function _inverse_vector(tab::AbstractMatrix{<:Integer})::Vector{Int}
  n = size(tab, 1)
  inv = zeros(Int, n)
  @inbounds for a in 1:n
    found = 0
    for b in 1:n
      if tab[a, b] == 1
        found = b
        break
      end
    end
    found == 0 && error("Group table has no inverse for element $a (no b with tab[a,b]=1).")
    inv[a] = found
  end
  return inv
end

"""
    TY_fusion_ring(tab; names=String[]) -> FusionRing

Build the Tambara–Yamagami fusion ring for a group with multiplication table `tab`.
Rank is n+1 (group elements + one extra object).
"""

function TY_fusion_ring(g::Group)::FusionRing
  return TY_fusion_ring(cayley_table(g))
end

function TY_fusion_ring(tab::AbstractMatrix{<:Integer}; names::Vector{String} = String[])
  _is_group_table(tab) || throw(
    ArgumentError(
      "FusionRingTY: tab must be a group multiplication table (identity=1, associative, latin square).",
    ),
  )
  n = size(tab, 1)
  r = n + 1

  mats = Matrix{Int}[]

  # For each simple object i=1..r, build its fusion matrix N_i.
  # This mirrors the Mathematica Which[...] table.
  @inbounds for i in 1:r
    Ni = zeros(Int, r, r)
    for j in 1:r
      if i <= n && j <= n
        k = tab[i, j]
        Ni[j, k] += 1
      elseif i <= n && j > n
        # group element ⊗ m = m
        Ni[j, r] += 1
      elseif i > n && j <= n
        # m ⊗ group element = m
        Ni[j, r] += 1
      else
        # m ⊗ m = sum_{g in G} g
        for k in 1:n
          Ni[j, k] += 1
        end
      end
    end
    push!(mats, Ni)
  end

  mt = _mats_to_mt(mats)

  # labels: 1..n are group elements, last is "m"
  labels = [string(i) for i in 1:n]
  push!(labels, "m")

  default_names = isempty(names) ? String[] : names
  return fusion_ring(mt; names = default_names, labels = labels)
end

# Helper function for TY_fusion_ring, HI_fusion_ring, and group_fusion_ring
"""
    _is_group_table(tab) -> Bool

Explicit check that `tab` is a group multiplication table on {1..n}
with identity element 1.
"""
function _is_group_table(tab::AbstractMatrix{<:Integer})::Bool
  n = size(tab, 1)
  size(tab, 2) == n || return false
  n ≥ 1 || return false

  # Entries in 1..n
  @inbounds for i in 1:n, j in 1:n
    x = tab[i, j]
    (1 <= x <= n) || return false
  end

  # Identity is 1
  @inbounds for i in 1:n
    tab[1, i] == i || return false
    tab[i, 1] == i || return false
  end

  # Latin square: each row/col is a permutation of 1..n
  seen = falses(n)
  @inbounds for i in 1:n
    fill!(seen, false)
    for j in 1:n
      seen[tab[i, j]] = true
    end
    all(seen) || return false

    fill!(seen, false)
    for j in 1:n
      seen[tab[j, i]] = true
    end
    all(seen) || return false
  end

  # Associativity
  @inbounds for i in 1:n, j in 1:n, k in 1:n
    tab[tab[i, j], k] == tab[i, tab[j, k]] || return false
  end

  return true
end

function cayley_table(grp::Group)::Matrix{Int64}
  els = collect(grp)
  return [findfirst(x -> x == g1 * g2, els) for g1 in els, g2 in els]
end

###############################################################################################
# Fusion Rings related to Quantum Groups
###############################################################################################

# TODO: add missing information
# PSU(2)_k
function psu2k_fusion_ring(k::Int)::FusionRing
  rk = div(k, 2) + 1
  mt = fill(0, rk, rk, rk)

  for a in 0:2:k, b in 0:2:k, c in 0:2:k
    if c ∈ range_psu2k(a, b, k)
      mt[div(a, 2) + 1, div(b, 2) + 1, div(c, 2) + 1] = 1
    else
      continue
    end
  end

  elnames = [
    "["*string(numerator((i-1)//2))*"/"*string(denominator((i-1)//2))*"]" for i in 1:rk
  ]

  return fusion_ring(mt; names = ["PSU(2)" * subscript_integer(k)], labels = elnames)
end

function range_psu2k(i::Int, j::Int, k::Int)
  return abs(i - j):2:min(i + j, 2k - i - j)
end

# TODO: add missing information
# SU(2)_k
function su2k_fusion_ring(k::Int)::FusionRing
  rk = k + 1
  mt = fill(0, rk, rk, rk)
  for a in 0:k, b in 0:k, c in 0:k
    if c ∈ range_psu2k(a, b, k)
      mt[a + 1, b + 1, c + 1] = 1
    else
      continue
    end
  end
  return fusion_ring(mt; names = ["SU(2)" * subscript_integer(k)], labels = string.(0:k))
end

# SO(2)_n fusion rings 

"""
    son2_fusion_ring(N::Int)::FusionRing

Return fusion ring ``\\text{SO}(N)_2`` (metaplectic) .
"""

#- odd `N`: uses `_son2_rules_odd(m)`
#- even `N ≡ 0 (mod 4)`: uses `_son2_rules_div4(N÷2)`
#- even `N ≡ 2 (mod 4)`: uses `_son2_rules_div2(N÷2)`

function son2_fusion_ring(N::Int)::FusionRing
  N ≥ 4 || throw(ArgumentError("son2_fusion_ring(N): requires integer N ≥ 4, got N=$N"))

  if isodd(N)
    mt     = _son2_rules_odd(N)
    labels = _son2_labels_odd(N)
  else
    p = N ÷ 2
    if N % 4 == 0
      mt = _son2_rules_div4(p)
    else
      mt = _son2_rules_div2(p)
    end
    labels = _son2_labels_even(p)
  end

  #  label count must match rank
  size(mt, 1) == length(labels) ||
    error("son2_fusion_ring: label length mismatch with mt rank")

  return fusion_ring(
    mt; names = ["SO($N)"*subscript_integer(2), "Metaplectic($N)"], labels = labels
  )
end

export metaplectic_fusion_ring

"""
    metaplectic_fusion_ring(N::Int)::FusionRing

Return fusion ring ``\\text{SO}(N)_2`` (metaplectic).
"""
metaplectic_fusion_ring(n::Int)::FusionRing = son2_fusion_ring(n)
metaplectic_fusion_ring(n::Int)::FusionRing = son2_fusion_ring(n)

function _son2_rules_odd(m::Integer)::Array{Int, 3}
  isodd(m) || throw(ArgumentError("_son2_rules_odd expects odd N, got N=$m"))
  m ≥ 5 || throw(ArgumentError("_son2_rules_odd expects N≥5 (odd), got N=$m"))

  r    = (m - 1) ÷ 2
  rank = (m + 7) ÷ 2

  # convenience
  ar(i) = _e(i, rank)

  mat1 = intidmat(rank)

  matZ = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 1
      ar(2)
    elseif i == 2
      ar(1)
    elseif 3 <= i <= 4
      ar(3 + mod(i, 2))   # i=3 -> 4, i=4 -> 3
    else
      ar(i)
    end
    matZ[i, :] .= v
  end

  # matXe1 = Table[ Which[...], {i,rank} ]
  matXe1 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 1
      ar(3)
    elseif i == 2
      ar(4)
    elseif i == 3
      ar(1) + sum(ar(j) for j in 5:rank)
    elseif i == 4
      ar(2) + sum(ar(j) for j in 5:rank)
    else
      ar(3) + ar(4)
    end
    matXe1[i, :] .= v
  end

  # matXe2 = Table[ Which[...], {i,rank} ]
  matXe2 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 1
      ar(4)
    elseif i == 2
      ar(3)
    elseif i == 3
      ar(2) + sum(ar(j) for j in 5:rank)
    elseif i == 4
      ar(1) + sum(ar(j) for j in 5:rank)
    else
      ar(3) + ar(4)
    end
    matXe2[i, :] .= v
  end

  # matY[j_] := Table[ Which[...], {i,rank} ]
  function matY(j::Int)::Matrix{Int}
    M = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
      v = if i == 1
        ar(j + 4)
      elseif i == 2
        ar(j + 4)
      elseif i == 3
        ar(3) + ar(4)
      elseif i == 4
        ar(3) + ar(4)
      else
        ii = i - 4
        if ii == j
          # ar[1] + ar[2] + ar[ Min[2j, m-2j] + 4 ]
          t = min(2*j, m - 2*j) + 4
          ar(1) + ar(2) + ar(t)
        else
          # ar[ Abs[ii-j] + 4 ] + ar[ Min[ii+j, m-ii-j+4] + 4 ]
          t1 = abs(ii - j) + 4
          t2 = min(ii + j, m - ii - j) + 4
          ar(t1) + ar(t2)
        end
      end
      M[i, :] .= v
    end
    return M
  end

  #   Transpose /@ Join[{mat1, matZ, matXe1, matXe2}, matY /@ Range[r]]
  mats = Matrix{Int}[]
  push!(mats, mat1)
  push!(mats, matZ)
  push!(mats, matXe1)
  push!(mats, matXe2)
  for j in 1:r
    push!(mats, matY(j))
  end

  # Convert fusion matrices to multiplication table
  return _mats_to_mt(mats)
end

# index convention (matches  Evaluate[...] = IdentityMatrix[rank]):
@inline 𝟙 = 1
@inline Θ = 2
@inline 𝟙 = 1
@inline Θ = 2
@inline Φ₁ = 3
@inline Φ₂ = 4
@inline σ₁ = 5
@inline σ₂ = 6
@inline τ₁ = 7
@inline τ₂ = 8
@inline Φ(j::Int) = 8 + j

# rulesdiv2[p_]
function _son2_rules_div2(p::Integer)::Array{Int, 3}
  p ≥ 1 || throw(ArgumentError("_son2_rules_div2 expects p≥1, got p=$p"))
  rank = p + 7
  maxphi = rank - 8  # = p-1

  ar(i) = _e(i, rank)

  # sums: sumEvenΛs = Σ_{i=2,4,...,p-1} Φ[i], sumOddΛs = Σ_{i=1,3,...,p-1} Φ[i]
  sumEven = zeros(Int, rank)
  sumOdd  = zeros(Int, rank)
  for i in 1:maxphi
    (isodd(i) ? (sumOdd[Φ(i)] += 1) : (sumEven[Φ(i)] += 1))
  end

  mats = Matrix{Int}[]

  # matId = IdentityMatrix[rank]
  matId = intidmat(rank)
  push!(mats, matId)

  # matΘ
  matTh = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(Θ)
    elseif i == Θ
      ar(𝟙)
    elseif i == Φ₁
      ar(Φ₂)
    elseif i == Φ₂
      ar(Φ₁)
    elseif i == σ₁
      ar(τ₁)
    elseif i == σ₂
      ar(τ₂)
    elseif i == τ₁
      ar(σ₁)
    elseif i == τ₂
      ar(σ₂)
    else
      ar(i) # Φ-lambdas fixed
    end
    matTh[i, :] .= v
  end
  push!(mats, matTh)

  # matΦ1
  matPhi1 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(Φ₁)
    elseif i == Θ
      ar(Φ₂)
    elseif i == Φ₁
      ar(Θ)
    elseif i == Φ₂
      ar(𝟙)
    elseif i == σ₁
      ar(σ₂)
    elseif i == σ₂
      ar(τ₁)
    elseif i == τ₁
      ar(τ₂)
    elseif i == τ₂
      ar(σ₁)
    else
      # Φ[p - (i-8)]
      j = i - 8
      ar(Φ(p - j))
    end
    matPhi1[i, :] .= v
  end
  push!(mats, matPhi1)

  # matΦ2
  matPhi2 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(Φ₂)
    elseif i == Θ
      ar(Φ₁)
    elseif i == Φ₁
      ar(𝟙)
    elseif i == Φ₂
      ar(Θ)
    elseif i == σ₁
      ar(τ₂)
    elseif i == σ₂
      ar(σ₁)
    elseif i == τ₁
      ar(σ₂)
    elseif i == τ₂
      ar(τ₁)
    else
      j = i - 8
      ar(Φ(p - j))
    end
    matPhi2[i, :] .= v
  end
  push!(mats, matPhi2)

  # matσ1
  matSig1 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(σ₁)
    elseif i == Θ
      ar(τ₁)
    elseif i == Φ₁
      ar(σ₂)
    elseif i == Φ₂
      ar(τ₂)
    elseif i == σ₁
      ar(Φ₂) + sumOdd
    elseif i == σ₂
      ar(𝟙) + sumEven
    elseif i == τ₁
      ar(Φ₁) + sumOdd
    elseif i == τ₂
      ar(Θ) + sumEven
    else
      # If[OddQ[i], σ2+τ2, σ1+τ1]
      isodd(i) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
    end
    matSig1[i, :] .= v
  end
  push!(mats, matSig1)

  # matσ2
  matSig2 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(σ₂)
    elseif i == Θ
      ar(τ₂)
    elseif i == Φ₁
      ar(τ₁)
    elseif i == Φ₂
      ar(σ₁)
    elseif i == σ₁
      ar(𝟙) + sumEven
    elseif i == σ₂
      ar(Φ₁) + sumOdd
    elseif i == τ₁
      ar(Θ) + sumEven
    elseif i == τ₂
      ar(Φ₂) + sumOdd
    else
      # If[EvenQ[i], σ2+τ2, σ1+τ1]
      iseven(i) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
    end
    matSig2[i, :] .= v
  end
  push!(mats, matSig2)

  # matτ1
  matTau1 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(τ₁)
    elseif i == Θ
      ar(σ₁)
    elseif i == Φ₁
      ar(τ₂)
    elseif i == Φ₂
      ar(σ₂)
    elseif i == σ₁
      ar(Φ₁) + sumOdd
    elseif i == σ₂
      ar(Θ) + sumEven
    elseif i == τ₁
      ar(Φ₂) + sumOdd
    elseif i == τ₂
      ar(𝟙) + sumEven
    else
      isodd(i) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
    end
    matTau1[i, :] .= v
  end
  push!(mats, matTau1)

  # matτ2
  matTau2 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(τ₂)
    elseif i == Θ
      ar(σ₂)
    elseif i == Φ₁
      ar(σ₁)
    elseif i == Φ₂
      ar(τ₁)
    elseif i == σ₁
      ar(Θ) + sumEven
    elseif i == σ₂
      ar(Φ₂) + sumOdd
    elseif i == τ₁
      ar(𝟙) + sumEven
    elseif i == τ₂
      ar(Φ₁) + sumOdd
    else
      iseven(i) ? ar(σ₂) + ar(τ₂) : ar(σ₁) .+ ar(τ₁)
    end
    matTau2[i, :] .= v
  end
  push!(mats, matTau2)

  # matΦ[j] for j = 1..(rank-8)
  function matPhi(j::Int)::Matrix{Int}
    M = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
      v = if i == 𝟙
        ar(Φ(j))
      elseif i == Θ
        ar(Φ(j))
      elseif i == Φ₁
        ar(Φ(p - j))
      elseif i == Φ₂
        ar(Φ(p - j))
      elseif i == σ₁
        isodd(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      elseif i == σ₂
        iseven(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      elseif i == τ₁
        isodd(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      elseif i == τ₂
        iseven(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      else
        ii = i - 8
        if ii == j && (2*j < p)
          ar(𝟙) + ar(Θ) + ar(Φ(2 * j))
        elseif ii == j && (2*j > p)
          ar(𝟙) + ar(Θ) + ar(Φ(2 * (p - j)))
        elseif ii + j < p
          ar(Φ(abs(ii - j))) + ar(Φ(ii + j))
        elseif ii + j > p
          ar(Φ(abs(ii - j))) + ar(Φ(2*p - ii - j))
        else
          # ii == p - j
          ar(Φ₁) + ar(Φ₂) + ar(Φ(abs(p - 2*ii)))
        end
      end
      M[i, :] .= v
    end
    return M
  end

  for j in 1:maxphi
    push!(mats, matPhi(j))
  end

  return _mats_to_mt(mats)
end

# rulesdiv4[p_]
function _son2_rules_div4(p::Integer)::Array{Int, 3}
  p ≥ 1 || throw(ArgumentError("_son2_rules_div4 expects p≥1, got p=$p"))
  rank = p + 7
  maxphi = rank - 8  # = p-1

  ar(i) = _e(i, rank)

  sumEven = zeros(Int, rank)
  sumOdd  = zeros(Int, rank)
  for i in 1:maxphi
    isodd(i) ? sumOdd[Φ(i)] += 1 : sumEven[Φ(i)] += 1
  end

  mats = Matrix{Int}[]

  matId = intidmat(rank)
  push!(mats, matId)

  # matΘ
  matTh = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(Θ)
    elseif i == Θ
      ar(𝟙)
    elseif i == Φ₁
      ar(Φ₂)
    elseif i == Φ₂
      ar(Φ₁)
    elseif i == σ₁
      ar(τ₁)
    elseif i == σ₂
      ar(τ₂)
    elseif i == τ₁
      ar(σ₁)
    elseif i == τ₂
      ar(σ₂)
    else
      ar(i)
    end
    matTh[i, :] .= v
  end
  push!(mats, matTh)

  # matΦ1
  matPhi1 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(Φ₁)
    elseif i == Θ
      ar(Φ₂)
    elseif i == Φ₁
      ar(𝟙)
    elseif i == Φ₂
      ar(Θ)
    elseif i == σ₁
      ar(τ₁)
    elseif i == σ₂
      ar(σ₂)
    elseif i == τ₁
      ar(σ₁)
    elseif i == τ₂
      ar(τ₂)
    else
      j = i - 8
      ar(Φ(p - j))
    end
    matPhi1[i, :] .= v
  end
  push!(mats, matPhi1)

  # matΦ2
  matPhi2 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(Φ₂)
    elseif i == Θ
      ar(Φ₁)
    elseif i == Φ₁
      ar(Θ)
    elseif i == Φ₂
      ar(𝟙)
    elseif i == σ₁
      ar(σ₁)
    elseif i == σ₂
      ar(τ₂)
    elseif i == τ₁
      ar(τ₁)
    elseif i == τ₂
      ar(σ₂)
    else
      j = i - 8
      ar(Φ(p - j))
    end
    matPhi2[i, :] .= v
  end
  push!(mats, matPhi2)

  # matσ1
  matSig1 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(σ₁)
    elseif i == Θ
      ar(τ₁)
    elseif i == Φ₁
      ar(τ₁)
    elseif i == Φ₂
      ar(σ₁)
    elseif i == σ₁
      ar(𝟙) + ar(Φ₂) + sumEven
    elseif i == σ₂
      sumOdd
    elseif i == τ₁
      ar(Θ) + ar(Φ₁) + sumEven
    elseif i == τ₂
      sumOdd
    else
      isodd(i) ? ar(σ₂) .+ ar(τ₂) : ar(σ₁) .+ ar(τ₁)
    end
    matSig1[i, :] .= v
  end
  push!(mats, matSig1)

  # matσ2
  matSig2 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(σ₂)
    elseif i == Θ
      ar(τ₂)
    elseif i == Φ₁
      ar(σ₂)
    elseif i == Φ₂
      ar(τ₂)
    elseif i == σ₁
      sumOdd
    elseif i == σ₂
      ar(𝟙) + ar(Φ₁) + sumEven
    elseif i == τ₁
      sumOdd
    elseif i == τ₂
      ar(Θ) + ar(Φ₂) + sumEven
    else
      iseven(i) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
    end
    matSig2[i, :] .= v
  end
  push!(mats, matSig2)

  # matτ1
  matTau1 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(τ₁)
    elseif i == Θ
      ar(σ₁)
    elseif i == Φ₁
      ar(σ₁)
    elseif i == Φ₂
      ar(τ₁)
    elseif i == σ₁
      ar(Θ) + ar(Φ₁) + sumEven
    elseif i == σ₂
      sumOdd
    elseif i == τ₁
      ar(𝟙) + ar(Φ₂) + sumEven
    elseif i == τ₂
      sumOdd
    else
      isodd(i) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
    end
    matTau1[i, :] .= v
  end
  push!(mats, matTau1)

  # matτ2
  matTau2 = zeros(Int, rank, rank)
  @inbounds for i in 1:rank
    v = if i == 𝟙
      ar(τ₂)
    elseif i == Θ
      ar(σ₂)
    elseif i == Φ₁
      ar(τ₂)
    elseif i == Φ₂
      ar(σ₂)
    elseif i == σ₁
      sumOdd
    elseif i == σ₂
      ar(Θ) + ar(Φ₂) + sumEven
    elseif i == τ₁
      sumOdd
    elseif i == τ₂
      ar(𝟙) + ar(Φ₁) + sumEven
    else
      iseven(i) ? ar(σ₂) .+ ar(τ₂) : ar(σ₁) .+ ar(τ₁)
    end
    matTau2[i, :] .= v
  end
  push!(mats, matTau2)

  # matΦ[j]
  function matPhi(j::Int)::Matrix{Int}
    M = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
      v = if i == 𝟙
        ar(Φ(j))
      elseif i == Θ
        ar(Φ(j))
      elseif i == Φ₁
        ar(Φ(p - j))
      elseif i == Φ₂
        ar(Φ(p - j))
      elseif i == σ₁
        isodd(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      elseif i == σ₂
        iseven(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      elseif i == τ₁
        isodd(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      elseif i == τ₂
        iseven(j) ? ar(σ₂) + ar(τ₂) : ar(σ₁) + ar(τ₁)
      else
        ii = i - 8
        if ii == j && (2*j < p)
          ar(𝟙) + ar(Θ) + ar(Φ(2*j))
        elseif ii == j && (2*j == p)
          ar(𝟙) + ar(Θ) + ar(Φ₁) + ar(Φ₂)
        elseif ii == j && (2*j > p)
          ar(𝟙) + ar(Θ) + ar(Φ(2*(p - j)))
        elseif ii + j < p
          ar(Φ(abs(ii - j))) + ar(Φ(ii + j))
        elseif ii + j > p
          ar(Φ(abs(ii - j))) + ar(Φ(2*p - ii - j))
        else
          # ii == p - j
          ar(Φ₁) + ar(Φ₂) + ar(Φ(abs(p - 2*ii)))
        end
      end
      M[i, :] .= v
    end
    return M
  end

  for j in 1:maxphi
    push!(mats, transpose(matPhi(j)))
  end

  return _mats_to_mt(mats)
end

# odd m: rank = (m+7)/2, elements are [1, Z, X_e1, X_e2, Y_1, ..., Y_r], r=(m-1)/2
function _son2_labels_odd(m::Int)::Vector{String}
  r = (m - 1) ÷ 2
  labels = String[bold_integer(1), "Z", "Xₑ₁", "Xₑ₂"]
  for j in 1:r
    push!(labels, "Y"*subscript_integer(j))
  end
  return labels
end

# even m: rank = p+7 with p=m/2, elements are [Id, Θ, Φ1, Φ2, σ1, σ2, τ1, τ2, Φ_1..Φ_{p-1}]
function _son2_labels_even(p::Int)::Vector{String}
  labels = String[bold_integer(1), "Θ", "Φ₁", "Φ₂", "σ₁", "σ₂", "τ₁", "τ₂"]
  for j in 1:(p - 1)
    push!(labels, "Φ"*subscript_integer(j))
  end
  return labels
end

# basis vector e_i in ℤ^rank
@inline function _e(i::Int, rank::Int)::Vector{Int}
  v = zeros(Int, rank)
  v[i] = 1
  return v
end

# convert a list of fusion matrices mats[a][b,c] into mt[a,b,c]
"""
    _mats_to_mt(mats) -> mt

Given mats[a] = N_a (rank×rank), return mt[a,b,c] = (N_a)[b,c].
"""
function _mats_to_mt(mats::Vector{<:AbstractMatrix{<:Integer}})::Array{Int, 3}
  r = length(mats)
  r ≥ 1 || error("_mats_to_mt: empty list of matrices")
  mt = zeros(Int, r, r, r)
  @inbounds for a in 1:r
    A = mats[a]
    size(A, 1) == r && size(A, 2) == r ||
      error("_mats_to_mt: mat $a has wrong size $(size(A)) (expected $r×$r)")
    mt[a, :, :] .= A
  end
  return mt
end
