# test/test_fixtures.jl
#
# Independently specified inputs and exact expected results.
# This file defines data providers; it does not run tests.

using FusionRings: fusion_ring
using Oscar: QQBar

# Conventions:
# - mt[a,b,c] is the coefficient of basis object c in a*b.
# - Basis index 1 is the unit.
# - Every fixture call returns fresh mutable arrays.
# - fixture_ring(f) constructs a validated ring without stored invariant caches.
# - Character rows may appear in a different order in computed output.
#   Compare rows up to permutation while preserving multiplicity.
# - Formal codegrees are sorted and retain repeated values.
# - proper_subrings excludes the unit-only subring and the whole ring.
# - grading_components gives the universal-grading partition of the basis.
# - upper_central_ranks includes the original ring and stops at a fixed point.
# - For noncommutative S3, characters/formal_codegrees = nothing indicates
#   unsupported package computations, not mathematical absence.
# - No fixture declares whether the ring is categorifiable.


# ============================================================
# Shared fixture construction
# ============================================================

function _fixture_data(
  name::Symbol,
  mt,
  labels,
  expected;
  cayley_table = nothing,
)
  return (
    name = name,
    input = (
      mt = copy(mt),
      labels = copy(labels),
      cayley_table =
        cayley_table === nothing ? nothing : copy(cayley_table),
    ),
    expected = merge(
      (
        rank = size(mt, 1),
        mt = copy(mt),
        labels = copy(labels),
      ),
      expected,
    ),
  )
end

function fixture_ring(fixture)
  return fusion_ring(
    copy(fixture.input.mt);
    labels = copy(fixture.input.labels),
    names = [String(fixture.name)],
    tex_names = [String(fixture.name)],
  )
end

# Encode a literal Cayley table directly.
# Expected tensors do not use group_fusion_ring or catalogue lookup.
function _fixture_group_mt(table)
  r = size(table, 1)

  return [
    Int(table[a, b] == c)
    for a in 1:r, b in 1:r, c in 1:r
  ]
end


# ============================================================
# Rank-one unit ring
# ============================================================

function fixture_unit()
  table = reshape([1], 1, 1)
  mt = _fixture_group_mt(table)

  expected = (
    duals = [1],
    commutative = true,
    group_ring = true,
    integral = true,
    weakly_integral = true,
    fpdims = QQBar.([1]),
    fpdim = QQBar(1),
    characters = reshape([QQBar(1)], 1, 1),
    formal_codegrees = QQBar.([1]),
    proper_subrings = Vector{Int}[],
    automorphisms = [[1]],
    adjoint_subset = [1],
    grading_components = [[1]],
    upper_central_ranks = [1],
    nilpotent = true,
  )

  return _fixture_data(
    :unit, mt, ["1"], expected;
    cayley_table = table,
  )
end


# ============================================================
# Z2: basis [1, g], with g^2 = 1
# ============================================================

function fixture_z2()
  table = [
    1 2
    2 1
  ]

  mt = _fixture_group_mt(table)

  expected = (
    duals = [1, 2],
    commutative = true,
    group_ring = true,
    integral = true,
    weakly_integral = true,
    fpdims = QQBar.([1, 1]),
    fpdim = QQBar(2),
    characters = QQBar.([
      1  1
      1 -1
    ]),
    formal_codegrees = QQBar.([2, 2]),
    proper_subrings = Vector{Int}[],
    automorphisms = [[1, 2]],
    adjoint_subset = [1],
    grading_components = [[1], [2]],
    upper_central_ranks = [2, 1],
    nilpotent = true,
  )

  return _fixture_data(
    :z2, mt, ["1", "g"], expected;
    cayley_table = table,
  )
end


# ============================================================
# Z3: basis [1, g, g^2]
# Non-self-dual pair and complex characters
# ============================================================

function fixture_z3()
  table = [
    1 2 3
    2 3 1
    3 1 2
  ]

  mt = _fixture_group_mt(table)
  omega = (-QQBar(1) + sqrt(QQBar(-3))) / 2

  expected = (
    duals = [1, 3, 2],
    commutative = true,
    group_ring = true,
    integral = true,
    weakly_integral = true,
    fpdims = QQBar.([1, 1, 1]),
    fpdim = QQBar(3),
    characters = QQBar.([
      1 1       1
      1 omega   omega^2
      1 omega^2 omega
    ]),
    formal_codegrees = QQBar.([3, 3, 3]),
    proper_subrings = Vector{Int}[],
    automorphisms = [
      [1, 2, 3],
      [1, 3, 2],
    ],
    adjoint_subset = [1],
    grading_components = [[1], [2], [3]],
    upper_central_ranks = [3, 1],
    nilpotent = true,
  )

  return _fixture_data(
    :z3, mt, ["1", "g", "g^2"], expected;
    cayley_table = table,
  )
end


# ============================================================
# Z4: basis [1, g, g^2, g^3]
# Proper subring [1, g^2] and a non-self-dual pair
# ============================================================

function fixture_z4()
  table = [
    1 2 3 4
    2 3 4 1
    3 4 1 2
    4 1 2 3
  ]

  mt = _fixture_group_mt(table)
  ii = sqrt(QQBar(-1))

  expected = (
    duals = [1, 4, 3, 2],
    commutative = true,
    group_ring = true,
    integral = true,
    weakly_integral = true,
    fpdims = QQBar.([1, 1, 1, 1]),
    fpdim = QQBar(4),
    characters = QQBar.([
      1  1     1  1
      1  ii   -1  (-ii)
      1 -1     1 -1
      1  (-ii) -1 ii
    ]),
    formal_codegrees = QQBar.([4, 4, 4, 4]),
    proper_subrings = [[1, 3]],
    automorphisms = [
      [1, 2, 3, 4],
      [1, 4, 3, 2],
    ],
    adjoint_subset = [1],
    grading_components = [[1], [2], [3], [4]],
    upper_central_ranks = [4, 1],
    nilpotent = true,
  )

  return _fixture_data(
    :z4, mt, ["1", "g", "g^2", "g^3"], expected;
    cayley_table = table,
  )
end


# ============================================================
# Klein four group: basis [1, a, b, ab]
# Multiple proper subrings and six automorphisms
# ============================================================

function fixture_v4()
  table = [
    1 2 3 4
    2 1 4 3
    3 4 1 2
    4 3 2 1
  ]

  mt = _fixture_group_mt(table)

  expected = (
    duals = [1, 2, 3, 4],
    commutative = true,
    group_ring = true,
    integral = true,
    weakly_integral = true,
    fpdims = QQBar.([1, 1, 1, 1]),
    fpdim = QQBar(4),
    characters = QQBar.([
      1  1  1  1
      1  1 -1 -1
      1 -1  1 -1
      1 -1 -1  1
    ]),
    formal_codegrees = QQBar.([4, 4, 4, 4]),
    proper_subrings = [
      [1, 2],
      [1, 3],
      [1, 4],
    ],
    automorphisms = [
      [1, 2, 3, 4],
      [1, 2, 4, 3],
      [1, 3, 2, 4],
      [1, 3, 4, 2],
      [1, 4, 2, 3],
      [1, 4, 3, 2],
    ],
    adjoint_subset = [1],
    grading_components = [[1], [2], [3], [4]],
    upper_central_ranks = [4, 1],
    nilpotent = true,
  )

  return _fixture_data(
    :v4, mt, ["1", "a", "b", "ab"], expected;
    cayley_table = table,
  )
end


# ============================================================
# Fibonacci: basis [1, tau], with tau^2 = 1 + tau
# ============================================================

function fixture_fibonacci()
  mt = zeros(Int, 2, 2, 2)

  mt[1, :, :] = [
    1 0
    0 1
  ]

  mt[2, :, :] = [
    0 1
    1 1
  ]

  s5 = sqrt(QQBar(5))
  phi = (1 + s5) / 2

  expected = (
    duals = [1, 2],
    commutative = true,
    group_ring = false,
    integral = false,
    weakly_integral = false,
    fpdims = QQBar.([1, phi]),
    fpdim = (5 + s5) / 2,
    characters = QQBar.([
      1 phi
      1 (1 - phi)
    ]),
    formal_codegrees = QQBar.([
      (5 - s5) / 2,
      (5 + s5) / 2,
    ]),
    proper_subrings = Vector{Int}[],
    automorphisms = [[1, 2]],
    adjoint_subset = [1, 2],
    grading_components = [[1, 2]],
    upper_central_ranks = [2],
    nilpotent = false,
  )

  return _fixture_data(
    :fibonacci, mt, ["1", "tau"], expected,
  )
end


# ============================================================
# Ising: basis [1, sigma, psi]
#
# sigma^2 = 1 + psi
# psi^2 = 1
# sigma*psi = psi*sigma = sigma
# ============================================================

function fixture_ising()
  mt = zeros(Int, 3, 3, 3)

  mt[1, :, :] = [
    1 0 0
    0 1 0
    0 0 1
  ]

  mt[2, :, :] = [
    0 1 0
    1 0 1
    0 1 0
  ]

  mt[3, :, :] = [
    0 0 1
    0 1 0
    1 0 0
  ]

  s2 = sqrt(QQBar(2))

  expected = (
    duals = [1, 2, 3],
    commutative = true,
    group_ring = false,
    integral = false,
    weakly_integral = true,
    fpdims = QQBar.([1, s2, 1]),
    fpdim = QQBar(4),
    characters = QQBar.([
      1 s2    1
      1 (-s2) 1
      1 0    -1
    ]),
    formal_codegrees = QQBar.([2, 4, 4]),
    proper_subrings = [[1, 3]],
    automorphisms = [[1, 2, 3]],
    adjoint_subset = [1, 3],
    grading_components = [[1, 3], [2]],
    upper_central_ranks = [3, 2, 1],
    nilpotent = true,
  )

  return _fixture_data(
    :ising, mt, ["1", "sigma", "psi"], expected,
  )
end


# ============================================================
# Representation ring of S3: basis [1, epsilon, rho]
#
# epsilon^2 = 1
# epsilon*rho = rho*epsilon = rho
# rho^2 = 1 + epsilon + rho
# ============================================================

function fixture_rep_s3()
  mt = zeros(Int, 3, 3, 3)

  mt[1, :, :] = [
    1 0 0
    0 1 0
    0 0 1
  ]

  mt[2, :, :] = [
    0 1 0
    1 0 0
    0 0 1
  ]

  mt[3, :, :] = [
    0 0 1
    0 0 1
    1 1 1
  ]

  expected = (
    duals = [1, 2, 3],
    commutative = true,
    group_ring = false,
    integral = true,
    weakly_integral = true,
    fpdims = QQBar.([1, 1, 2]),
    fpdim = QQBar(6),
    characters = QQBar.([
      1  1  2
      1 -1  0
      1  1 -1
    ]),
    formal_codegrees = QQBar.([2, 3, 6]),
    proper_subrings = [[1, 2]],
    automorphisms = [[1, 2, 3]],
    adjoint_subset = [1, 2, 3],
    grading_components = [[1, 2, 3]],
    upper_central_ranks = [3],
    nilpotent = false,
  )

  return _fixture_data(
    :rep_s3, mt, ["1", "epsilon", "rho"], expected,
  )
end


# ============================================================
# Noncommutative S3 group ring
#
# Basis [1, r, r^2, s, rs, r^2s]
# r^3 = s^2 = 1 and s*r = r^(-1)*s
# ============================================================

function fixture_s3()
  table = [
    1 2 3 4 5 6
    2 3 1 5 6 4
    3 1 2 6 4 5
    4 6 5 1 3 2
    5 4 6 2 1 3
    6 5 4 3 2 1
  ]

  mt = _fixture_group_mt(table)

  expected = (
    duals = [1, 3, 2, 4, 5, 6],
    commutative = false,
    group_ring = true,
    integral = true,
    weakly_integral = true,
    fpdims = QQBar.([1, 1, 1, 1, 1, 1]),
    fpdim = QQBar(6),
    characters = nothing,
    formal_codegrees = nothing,
    proper_subrings = [
      [1, 2, 3],
      [1, 4],
      [1, 5],
      [1, 6],
    ],
    automorphisms = [
      [1, 2, 3, 4, 5, 6],
      [1, 2, 3, 5, 6, 4],
      [1, 2, 3, 6, 4, 5],
      [1, 3, 2, 4, 6, 5],
      [1, 3, 2, 5, 4, 6],
      [1, 3, 2, 6, 5, 4],
    ],
    adjoint_subset = [1],
    grading_components = [[1], [2], [3], [4], [5], [6]],
    upper_central_ranks = [6, 1],
    nilpotent = true,
  )

  return _fixture_data(
    :s3,
    mt,
    ["1", "r", "r^2", "s", "rs", "r^2s"],
    expected;
    cayley_table = table,
  )
end


# ============================================================
# Rank-two multiplicity-two ring: basis [1, x]
#
# x^2 = 1 + 2x
# This fixture makes no categorifiability claim.
# ============================================================

function fixture_rank2_m2()
  mt = zeros(Int, 2, 2, 2)

  mt[1, :, :] = [
    1 0
    0 1
  ]

  mt[2, :, :] = [
    0 1
    1 2
  ]

  s2 = sqrt(QQBar(2))

  expected = (
    duals = [1, 2],
    commutative = true,
    group_ring = false,
    integral = false,
    weakly_integral = false,
    fpdims = QQBar.([1, 1 + s2]),
    fpdim = 4 + 2*s2,
    characters = QQBar.([
      1 (1 + s2)
      1 (1 - s2)
    ]),
    formal_codegrees = QQBar.([
      4 - 2*s2,
      4 + 2*s2,
    ]),
    proper_subrings = Vector{Int}[],
    automorphisms = [[1, 2]],
    adjoint_subset = [1, 2],
    grading_components = [[1, 2]],
    upper_central_ranks = [2],
    nilpotent = false,
  )

  return _fixture_data(
    :rank2_m2, mt, ["1", "x"], expected,
  )
end


# ============================================================
# Fixture collection
# ============================================================

"""Return all ten fixtures, each with fresh arrays and exact expectations."""
function small_ring_fixtures()
  return (
    fixture_unit(),
    fixture_z2(),
    fixture_z3(),
    fixture_z4(),
    fixture_v4(),
    fixture_fibonacci(),
    fixture_ising(),
    fixture_rep_s3(),
    fixture_s3(),
    fixture_rank2_m2(),
  )
end


# ============================================================
# Invalid constructor inputs
# ============================================================

"""
    invalid_ring_inputs()

Return raw inputs for focused constructor-rejection tests.

Each case provides:
- name
- mt
- kwargs
- exception
- message_contains

No invalid FusionRing is constructed here.

The nonassociative example has nonnegative integer entries, a two-sided unit,
and unique duals, so it reaches the associativity check.
"""
function invalid_ring_inputs()
  z2 = fixture_z2().input.mt

  negative = copy(z2)
  negative[2, 2, 1] = -1

  missing_unit = copy(z2)
  missing_unit[1, 2, 2] = 0

  missing_dual = copy(z2)
  missing_dual[2, 2, 1] = 0
  missing_dual[2, 2, 2] = 1

  # Basis [1, x, y], with:
  # x^2 = 1+y
  # xy = yx = x+y
  # y^2 = 1+x
  #
  # The coefficient of y in (x*x)*y is 1.
  # The coefficient of y in x*(x*y) is 2.
  nonassociative = zeros(Int, 3, 3, 3)

  nonassociative[1, :, :] = [
    1 0 0
    0 1 0
    0 0 1
  ]

  nonassociative[2, :, :] = [
    0 1 0
    1 0 1
    0 1 1
  ]

  nonassociative[3, :, :] = [
    0 0 1
    0 1 1
    1 1 0
  ]

  return (
    (
      name = :rank_zero,
      mt = zeros(Int, 0, 0, 0),
      kwargs = (;),
      exception = ErrorException,
      message_contains = "equal side lengths",
    ),
    (
      name = :noncubic,
      mt = zeros(Int, 2, 2, 3),
      kwargs = (;),
      exception = ErrorException,
      message_contains = "equal side lengths",
    ),
    (
      name = :negative_coefficient,
      mt = negative,
      kwargs = (;),
      exception = ErrorException,
      message_contains = "non-negative integers",
    ),
    (
      name = :missing_unit,
      mt = missing_unit,
      kwargs = (;),
      exception = ErrorException,
      message_contains = "act as unit object",
    ),
    (
      name = :missing_dual,
      mt = missing_dual,
      kwargs = (;),
      exception = ErrorException,
      message_contains = "unique inverse",
    ),
    (
      name = :nonassociative,
      mt = nonassociative,
      kwargs = (;),
      exception = ErrorException,
      message_contains = "associativity",
    ),
    (
      name = :wrong_label_count,
      mt = copy(z2),
      kwargs = (labels = ["1"],),
      exception = ErrorException,
      message_contains = "labels length",
    ),
    (
      # A valid mathematical table stored in an unsupported element type.
      name = :float_storage,
      mt = Float64.(z2),
      kwargs = (;),
      exception = MethodError,
      message_contains = nothing,
    ),
  )
end
