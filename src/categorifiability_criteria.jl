export csp_criterion,
  pdc_criterion,
  dn_criterion,
  is_d_number,
  ecc_criterion,
  lagrange_criterion,
  zsc_criterion,
  osc_criterion

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                     commutative schur product criterion                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#=
@article{LIU2021107905,
title = {Fusion bialgebras and Fourier analysis: Analytic obstructions for unitary categorification},
journal = {Advances in Mathematics},
volume = {390},
pages = {107905},
year = {2021},
issn = {0001-8708},
doi = {https://doi.org/10.1016/j.aim.2021.107905},
url = {https://www.sciencedirect.com/science/article/pii/S0001870821003443},
author = {Zhengwei Liu and Sebastien Palcoux and Jinsong Wu},
keywords = {Quantum Fourier analysis, Subfactors, Planar algebras, Fusion rings, Unitary categorification},
abstract = {We introduce fusion bialgebras and their duals and systematically study their Fourier analysis. As an application, we discover new efficient analytic obstructions on the unitary categorification of fusion rings. We prove the Hausdorff-Young inequality, uncertainty principles for fusion bialgebras and their duals. We show that the Schur product property, Young's inequality and the sum-set estimate hold for fusion bialgebras, but not always on their duals. If the fusion ring is the Grothendieck ring of a unitary fusion category, then these inequalities hold on the duals. Therefore, these inequalities are analytic obstructions of categorification. We classify simple integral fusion rings of Frobenius type up to rank 8 and of Frobenius-Perron dimension less than 4080. We find 34 ones, 4 of which are group-like and 28 of which can be eliminated by applying the Schur product property on the dual. In general, these inequalities are obstructions to subfactorize fusion bialgebras.}
}

Corollary 8.5
=#

"""
csp_criterion(ring; force_compute=false)

Return `true` when the commutative Schur product criterion obstructs a
unitary categorification. A noncommutative ring returns `false` because this
commutative criterion does not apply.
"""

# characters(ring) stores characters by row and simple objects by column, so
# the theorem's `λ[i,j]` is represented by `chars[j,i]` here.
 
  
function csp_criterion(
  ring::FusionRing;
  force_compute::Bool = false,
)::Bool
  return first(
    csp_criterion_explanation(
      ring;
      force_compute = force_compute,
    ),
  )
end





#changed: Always return   explanation tuple when 
# criterion does not apply or finds no obstruction.
function csp_criterion_explanation(
  ring::FusionRing;
  force_compute::Bool = false,
)::Tuple{
  Bool,
  Union{Nothing,Tuple{QQBarFieldElem,NTuple{3,Int}}},
  String
}
  if !is_commutative(ring)
    return (
      false,
      nothing,
      "The commutative Schur product criterion does not apply to noncommutative rings.",
    )
  end

  chars = characters(ring; force_compute = force_compute)
  dimensions = fpdims(ring; force_compute = force_compute)
  r = rank(ring)

  for j1 in 1:r, j2 in 1:r, j3 in 1:r
    coeff = sum(
      chars[j1, i] *
      chars[j2, i] *
      chars[j3, i] /
      dimensions[i]
      for i in 1:r
    )

    if !is_real(coeff) || coeff < 0
      indices = (j1, j2, j3)

      return (
        true,
        (coeff, indices),
        csp_explanation(coeff, indices),
      )
    end
  end

  return (
    false,
    nothing,
    "The commutative Schur product criterion found no obstruction.",
  )
end


#changed: Format  exact | approximate coefficient correctly and remove
#  undefined aprstr variable.
function csp_explanation(
  coeff::QQBarFieldElem,
  indices::NTuple{3,Int},
)::String
  j1, j2, j3 = indices

  exact_string = qqb_id_rep(qqb_id(coeff))
  approximate_string = _format_complex_float(
    string(ComplexF64(coeff)),
  )

  return "The fusion ring does not have a unitary categorification " *
         "because \$s = \\sum_i " *
         "\\frac{" *
         "\\lambda_{$(j1),i}" *
         "\\lambda_{$(j2),i}" *
         "\\lambda_{$(j3),i}" *
         "}{\\lambda_{1,i}}\$ = " *
         exact_string *
         " \\approx " *
         approximate_string *
         ", but \$s \\geq 0\$ according to the commutative Schur product " *
         "criterion from 10.1016/j.aim.2021.107905, Corollary 8.5."
end


function _format_complex_float(s::String)::String
  s = replace(s, " + 0.0im" => "")
  return replace(s, "im" => "i")
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                       pivotal drinfeld center criterion                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# pdc_criterion returns true if ring has no complex pivotal categorification due to the pivotal Drinfeld center criterion

"""
pdc_criterion(fr; force_compute=false)

Return `true` when the pivotal Drinfeld-center criterion obstructs a complex
pivotal categorification.
"""
#changed: Implement the formal-codegree ratio test with exact integrality, force_compute, and obstruction semantics.
function pdc_criterion(
  fr::FusionRing; force_compute::Bool = false
)::Bool
  is_commutative(fr) || return false

  codegrees = formal_codegrees(fr; force_compute = force_compute)
  for candidate in codegrees
    all(codegree -> is_algebraic_integer(candidate/codegree), codegrees) &&
      return false
  end

  return true
end



#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                    pseudo-unitary drinfeld center criterion                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

function pudc_criterion(
  fr::FusionRing; force_compute::Bool = false
)::Bool
  first( pdc_criterion_explanation(fr,force_compute=force_compute) )
end


function pudc_criterion_explanation(
  fr::FusionRing; force_compute::Bool = false
)::Tuple{Bool,Tuple{QQBarFieldElem,Int64},String}
  is_commutative(fr) || return false

  codegrees = formal_codegrees(fr; force_compute = force_compute)

  for i in 1:rank(fr)
    frac = codegrees[1]//codegrees[i]
    if !is_algebraic_integer(frac)
      return ( true, (frac, i), pdc_explanation(frac,i) )
    end
  end

  return false
end

  
function pudc_explanation( frac, int )
  qqbstr = qqb_id(frac)
  
  "The fusion ring does not have a pseudo-untary categorification"*
  " because \$\\frac{c_1}{c_$(int)}\$ ="*qqbstr*
  ", where c_i is the i'th formal codegree, is not an algebraic"*
  " integer which it needs to be according to the pseudo-unitary"*
  " version of the Drinfeld center criterion from 10.1007/s11005-022-01542-1 "*
  "Theorem 2.4."
end
  
#function pudc_criterion(fr::FusionRing)
#  !is_commutative(fr) && return false
#
#  chars = characters(fr)
#
#end
#  Returns True if ring has no complex pseudo-unitary categorification
#
#
#PackageExport["PUDCC"]
#
#PUDCC[ ring_FusionRing?CommutativeQ ] :=
#  With[{ chars = FusionRingCharacters[ring] },
#    And @@ Flatten @
#    AlgebraicIntegerQ[ chars[[1]].ConjugateTranspose[chars[[1]]] / chars ]
#  ];
#
#
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                               D-number criterion                                ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
# Source:https://arxiv.org/pdf/0810.3242
#
# The function returns true if the fusion ring has no complex categorification
#
"""
dn_criterion(fr; force_compute=false)

Return `true` when at least one formal codegree is not a d-number, obstructing
a complex categorification. This executable test currently returns `false` for
noncommutative rings because the package's general noncommutative formal-codegree
path does not yet recover the irreducible-representation dimensions needed to
separate `f_E` from `f_E * dim(E)`.
"""
#changed: Apply Ostrik's d-number theorem to every supported commutative formal codegree; any failure obstructs categorification.
function dn_criterion(
  fr::FusionRing; force_compute::Bool = false
)::Bool
  is_commutative(fr) || return false
  codegrees = formal_codegrees(fr; force_compute = force_compute)
  return any(codegree -> !is_d_number(codegree), codegrees)
end

"""Return whether the algebraic number `z` is a d-number."""
#changed: Require algebraic integrality and correct the minimal-polynomial coefficient/divisibility direction.
function is_d_number(z::QQBarFieldElem)::Bool
  is_algebraic_integer(z) || return false

  polynomial_ring_qq, _ = polynomial_ring(QQ, :x)
  polynomial = minpoly(polynomial_ring_qq, z)
  n = degree(polynomial)
  n <= 1 && return true

  constant_term = coeff(polynomial, 0)
  constant_term == 0 && return false

  for i in 1:(n - 1)
    coefficient = coeff(polynomial, n - i)
    denominator(coefficient^n / constant_term^i) == 1 || return false
  end
  return true
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                          extended cyclotomic criterion                          ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""
ecc_criterion(fr::FusionRing)

Return true if  extended cyclotomic criterion obstructs a complex
categorification of `fr`.

For every fusion matrix,  function computes its minimal polynomial over
the rationals and then the Galois group of its splitting field. A non-abelian
Galois group obstructs complex categorification.

This currently only works for commutative rings, a non-commutative ring will automatically return false.
"""

function _has_nonabelian_splitting_field(polynomial)::Bool
  # Constant and linear polynomials already split over QQ.
  degree(polynomial) <= 1 && return false

  galois_group_of_polynomial, _ = galois_group(polynomial)
  return !is_abelian(galois_group_of_polynomial)
end

#changed: Implement the extended cyclotomic obstruction using the exact
# minimal polynomial and Galois group of each fusion matrix.
function ecc_criterion(fr::FusionRing)::Bool
  is_commutative(fr) || return false

  polynomial_ring_qq, _ = polynomial_ring(QQ, :x)
  multiplication = multiplication_table(fr)

  for i in 1:rank(fr)
    fusion_matrix = matrix(QQ, multiplication[i, :, :])
    polynomial = minpoly(polynomial_ring_qq, fusion_matrix)

    _has_nonabelian_splitting_field(polynomial) && return true
  end

  return false
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                Lagrange criterion                               ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""
lagrange_criterion(fr::FusionRing; force_compute=false)

Returns true if  Lagrange criterion obstructs categorification.

For every proper, nontrivial fusion-closed based subring `S`, the function
checks whether

    FPdim(fr) / FPdim(S)

is an algebraic integer. If any quotient is not an algebraic integer, the
fusion ring cannot be categorifiable
"""

#changed: Implement the Lagrange obstruction using fusion-closed subsets and
# the parent ring's FP dimensions, without constructing separate FusionRings.
function lagrange_criterion(
  fr::FusionRing; force_compute::Bool = false
)::Bool
  dimensions = fpdims(fr; force_compute = force_compute)
  squared_dimensions = dimensions .^ 2
  total_dimension = sum(squared_dimensions)

  for subset in sub_fusion_ring_subsets(fr)
    subring_dimension = sum(i -> squared_dimensions[i], subset)
    quotient = total_dimension / subring_dimension

    !is_algebraic_integer(quotient) && return true
  end

  return false
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                      spectrum criterion helper functions                        ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

#changed: Build indexed candidate lists so ZSC  OSC do not repeatedly scan
# every nonzero structure constant for fixed-coordinate matches.
function nonzero_structure_constant_lookups(nonzero, r::Int)
  tuple_type = eltype(nonzero)

  by_second = [tuple_type[] for _i in 1:r]
  by_second_third = [
    tuple_type[] for _i in 1:r, _j in 1:r
  ]
  by_third = [tuple_type[] for _i in 1:r]
  by_first_second = [
    tuple_type[] for _i in 1:r, _j in 1:r
  ]

  for ind in nonzero
    first_index, second_index, third_index = ind

    push!(by_second[second_index], ind)
    push!(by_second_third[second_index, third_index], ind)
    push!(by_third[third_index], ind)
    push!(by_first_second[first_index, second_index], ind)
  end

  return (
    by_second,
    by_second_third,
    by_third,
    by_first_second,
  )
end


#changed: Test  three crit1 sums one at a time and stop any sum as soon as
# it exceeds one.
function crit1_has_one(i, d, mt)::Bool
  length(i) == 6 || throw(ArgumentError("crit1 requires six indices"))

  r = size(mt, 1)

  total = 0
  for k in 1:r
    total += mt[i[5], i[4], k] * mt[i[3], d[i[1]], k]
    total > 1 && break
  end
  total == 1 && return true

  total = 0
  for k in 1:r
    total += mt[i[2], d[i[4]], k] * mt[i[3], d[i[6]], k]
    total > 1 && break
  end
  total == 1 && return true

  total = 0
  for k in 1:r
    total += mt[d[i[5]], i[2], k] * mt[i[6], d[i[1]], k]
    total > 1 && return false
  end

  return total == 1
end


#changed: same as before
function crit2_has_one(i, d, mt)::Bool
  length(i) == 6 || throw(ArgumentError("crit2 requires six indices"))

  r = size(mt, 1)

  total = 0
  for k in 1:r
    total += mt[i[2], i[4], k] * mt[i[3], d[i[6]], k]
    total > 1 && break
  end
  total == 1 && return true

  total = 0
  for k in 1:r
    total += mt[i[5], d[i[4]], k] * mt[i[3], d[i[1]], k]
    total > 1 && break
  end
  total == 1 && return true

  total = 0
  for k in 1:r
    total += mt[d[i[2]], i[5], k] * mt[i[1], d[i[6]], k]
    total > 1 && return false
  end

  return total == 1
end


#changed:  ZSC stops as soon as  crit3 sum  known to be nonzero.
function crit3_is_zero(i, d, mt)::Bool
  length(i) == 6 || throw(ArgumentError("crit3 requires six indices"))

  r = size(mt, 1)

  for k in 1:r
    product =
      mt[i[1], i[4], k] *
      mt[d[i[2]], i[5], k] *
      mt[i[3], d[i[6]], k]

    product != 0 && return false
  end

  return true
end


#changed:  OSC stops as soon as  nonnegative crit3 sum > 1
function crit3_is_one(i, d, mt)::Bool
  length(i) == 6 || throw(ArgumentError("crit3 requires six indices"))

  r = size(mt, 1)
  total = 0

  for k in 1:r
    total +=
      mt[i[1], i[4], k] *
      mt[d[i[2]], i[5], k] *
      mt[i[3], d[i[6]], k]

    total > 1 && return false
  end

  return total == 1
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              zero spectrum criterion                            ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""
zsc_criterion(ring)

Returns true if the Zero Spectrum Criterion rules out categorifiability.
"""
#changed: Search only indexed matching structure constants, begin from entries
# equal to one, and use short-circuiting.
function zsc_criterion(ring::FusionRing)::Bool
  mt = multiplication_table(ring)
  r = size(mt, 1)

  d = [conjugate_element(ring, i) for i in 1:r]
  nonzero = nonzero_structure_constants(ring)

  (
    by_second,
    by_second_third,
    by_third,
    by_first_second,
  ) = nonzero_structure_constant_lookups(nonzero, r)

  # Each position is a CartesianIndex(i2, i1, i3) satisfying
  # mt[i2, i1, i3] == 1.
  for position in findall(==(1), mt)
    i2, i1, i3 = Tuple(position)

    # Tuples of the form (i4, i1, i6).
    for ind1 in by_second[i1]
      i4 = ind1[1]
      i6 = ind1[3]

      # Tuples of the form (i5, i4, i2).
      for ind2 in by_second_third[i4, i2]
        i5 = ind2[1]

        mt[i5, i6, i3] != 0 || continue

        crit1_has_one(
          (i1, i2, i3, i4, i5, i6),
          d,
          mt,
        ) || continue

        # Tuples of the form (i7, i9, i1).
        for ind3 in by_third[i1]
          i7 = ind3[1]
          i9 = ind3[2]

          # Tuples of the form (i2, i7, i8).
          for ind4 in by_first_second[i2, i7]
            i8 = ind4[3]

            mt[i8, i9, i3] != 0 || continue

            crit3_is_zero(
              (i4, i5, i6, i7, i8, i9),
              d,
              mt,
            ) || continue

            crit2_has_one(
              (i1, i2, i3, i7, i8, i9),
              d,
              mt,
            ) || continue

            return true
          end
        end
      end
    end
  end

  return false
end


#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              one spectrum criterion                             ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""
osc_criterion(ring)

Returns true if  One Spectrum Criterion rules out categorifiability.
"""
#changed: Search only indexed matching structure constants and use
# short-circuiting spectrum predicates.
function osc_criterion(ring::FusionRing)::Bool
  mt = multiplication_table(ring)
  r = size(mt, 1)

  d = [conjugate_element(ring, i) for i in 1:r]
  nonzero = nonzero_structure_constants(ring)

  (
    by_second,
    by_second_third,
    by_third,
    by_first_second,
  ) = nonzero_structure_constant_lookups(nonzero, r)

  # Zero entries are usually common, so retaining the direct traversal avoids
  # allocating a potentially large list containing every zero position.
  for i2 in 1:r, i1 in 1:r, i3 in 1:r
    mt[i2, i1, i3] == 0 || continue

    # Tuples of the form (i4, i1, i6).
    for ind1 in by_second[i1]
      i4 = ind1[1]
      i6 = ind1[3]

      # Tuples of the form (i5, i4, i2).
      for ind2 in by_second_third[i4, i2]
        i5 = ind2[1]

        mt[i5, i6, i3] != 0 || continue

        # Tuples of the form (i7, i9, i1).
        for ind3 in by_third[i1]
          i7 = ind3[1]
          i9 = ind3[2]

          for i0 in 1:r
            mt[i4, i7, i0] == 1 || continue
            mt[i6, d[i9], i0] == 1 || continue

            crit1_has_one(
              (i9, i0, i6, i7, i4, i1),
              d,
              mt,
            ) || continue

            # These tuples already guarantee mt[i2, i7, i8] != 0.
            for ind4 in by_first_second[i2, i7]
              i8 = ind4[3]

              mt[i8, i9, i3] != 0 || continue
              mt[d[i5], i8, i0] == 1 || continue

              crit3_is_one(
                (i4, i5, i6, i7, i8, i9),
                d,
                mt,
              ) || continue

              crit1_has_one(
                (i7, i2, i8, i4, i5, i0),
                d,
                mt,
              ) || continue

              crit1_has_one(
                (i9, i8, i3, i0, i5, i6),
                d,
                mt,
              ) || continue

              return true
            end
          end
        end
      end
    end
  end

  return false
end
