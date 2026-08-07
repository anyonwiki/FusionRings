#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                     commutative schur product criterion                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# TODO: check whether s is correct
#Todo
# csp_criterion( fusion_ring ) returns true if fusion_ring does not have a unitary categorification due to the commutative schur product criterion

function csp_criterion(ring::FusionRing)
  is_commutative(ring) || return false

  chars = characters(ring)
  r = rank(ring)

  for j1 in 1:r, j2 in 1:r, j3 in 1:r
    s = sum(chars[i, j1] * chars[i, j2] * chars[i, j3] / chars[i, 1] for i in 1:r)
    if is_real(s) && s < 0
      return true
    else
      continue
    end
  end

  return false
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                       pivotal drinfeld center criterion                         ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

# pdc_criterion returns true if ring has no complex pivotal categorification due to the pivotal Drinfeld center criterion

function pdc_criterion(fr::FusionRing)::Bool
  !is_commutative(fr) && return false

  r = rank(fr)
  fcds = formal_codegrees(fr)

  for j in 1:r
    ratios = [fcds[j]/fcds[i] for i in 1:r]
    if all(is_algebraic_integer, ratios)
      return false
    end
  end

  return true
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                    pseudo-unitary drinfeld center criterion                     ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

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
export dn_criterion

function dn_criterion(fr::FusionRing)::Bool
  !is_commutative(fr) && return false

  fcds = formal_codegrees(fr)

  for z in fcds
    !is_d_number(z) && return true
  end

  return false
end

function is_d_number(z::QQBarFieldElem)::Bool
  R, _ = polynomial_ring(QQ, :x)
  mp   = minpoly(R, z)
  n    = degree(mp)

  an  = constant_coefficient(mp)
  cfs = collect(coefficients(mp))[2:(end - 1)]

  is_empty(cfs) && return true

  divint(i) = is_integer((an^i) / (cfs[i]^n))

  return all(divint, 1:(n - 1))
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                          extended cyclotomic criterion                          ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#  The function returns True if the fusion ring cannot be categorified.
#
#
#
#
#
#
#
#
#
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                                Lagrange criterion                               ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#  The function returns False if the fusion ring cannot be categorified.
#
#
#PackageExport["LCriterion"]
#
#LCriterion[ ring_FusionRing ] :=
#  Module[{ subringDims },
#    subringDims =
#      FPDim /@
#      DeleteDuplicates @
#      SubFusionRings[ ring ][[;;,2]];
#
#    Not[And @@ AlgebraicIntegerQ @ RootReduce[ FPDim[ring]/subringDims ]]
#
#  ];
#
#
#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              zero spectrum criterion                            ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#  The function returns True if the fusion ring cannot be categorified.
#
#PackageExport["ZSCriterion"]
#
#ZSCriterion::usage =
#  "ZSCriterion[ fusionRing ] returns True if the fusion ring fusionRing cannot be categorified due to " <>
#  "the Zero Spectrum criterion. (arXiv:2203.06522v1)";
#
#(*See arXiv:2203.06522v1 for more info.";*)
#
#SetAttributes[ ZSC, Listable ];
#
#ZSCriterion[ ring_FusionRing ] :=
#  Module[{ mt, non0Cons, ones, d, i1, i2, i3, i4, i5, i6, i7, i8, i9, matches1, matches2, matches3, matches4 },
#    mt =
#      MT[ring];
#    non0Cons =
#      NonZeroStructureConstants[ring];
#    ones =
#      Position[ mt, x_Integer/; x == 1 ];
#    d =
#      CC[ring] /@ Range[Rank[ring]];
#
#    Catch[
#      Do[
#        { i2, i1, i3 } = ind;
#        matches1 = Cases[ non0Cons, { i4_, i1, i6_ } ];
#        Do[
#          { i4, i6 } = ind1[[{1,3}]];
#          matches2 = Cases[ non0Cons, { i5_, i4, i2 } ];
#          Do[
#            i5 = ind2[[1]];
#            If[
#              mt[[ i5, i6, i3 ]] != 0 &&
#              MemberQ[ crit1[ {i1,i2,i3,i4,i5,i6}, d, mt ], 1 ],
#              matches3 = Cases[ non0Cons, { i7_, i9_, i1 } ];
#              Do[
#                { i7, i9 } = ind3[[{1,2}]];
#                matches4 = Cases[ non0Cons, { i2, i7, i8_ } ];
#                Do[
#                  i8 = ind4[[3]];
#                  If[
#                    mt[[ i8, i9, i3 ]] != 0 &&
#                    crit3[ { i4, i5, i6, i7, i8, i9 }, d, mt ] == 0 &&
#                    MemberQ[ crit2[ { i1, i2, i3, i7, i8, i9 }, d, mt ], 1 ],
#                    Throw[ True ]
#                  ]
#                  ,{ ind4, matches4 }]
#                ,{ ind3, matches3 }]
#            ]
#            ,{ ind2, matches2 }]
#          ,{ ind1, matches1 }]
#        , { ind, ones } ];
#      False
#    ]
#
#  ];
#
function crit1(i::Vector{Int64}, d::Vector{Int64}, mt::Array{Int64, 3})::Bool
  r = size(mt, 1)
  return sum(mt[i[5], i[4], k] * mt[i[3], d[i[1]], k] for k in 1:r) == 1 ||
         sum(mt[i[2], d[i[4]], k] * mt[i[3], d[i[6]], k] for k in 1:r) == 1 ||
         sum(mt[d[i[5]], i[2], k] * mt[i[6], d[i[1]], k] for k in 1:r) == 1
end

function crit2(i::Vector{Int64}, d::Vector{Int64}, mt::Array{Int64, 3})::Bool
  r = size(mt, 1)
  return sum(mt[i[2], i[4], k] * mt[i[3], d[i[6]], k] for k in 1:r) == 1 ||
         sum(mt[i[5], d[i[4]], k] * mt[i[3], d[i[1]], k] for k in 1:r) == 1 ||
         sum(mt[d[i[2]], i[5], k] * mt[i[1], d[i[6]], k] for k in 1:r) == 1
end

function crit3(i::Vector{Int64}, d::Vector{Int64}, mt::Array{Int64, 3})::Bool
  r = size(mt, 1)
  return sum(
    mt[i[1], i[4], k] * mt[d[i[2]], i[5], k] * mt[i[3], d[i[6]], k] for k in 1:r
  ) == 0
end

#┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#┃                              one spectrum criterion                             ┃
#┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#PackageExport["OSCriterion"]
#
#OSCriterion::usage =
#  "OSCriterion[ fusionRing ] returns True if the fusion ring fusionRing cannot be categorified due to " <>
#  "the One Spectrum criterion. (arXiv:2203.06522v1)";
#
#OSCriterion[ ring_FusionRing ] :=
#  Module[{ mt, non0Cons, zeros, d, i1, i2, i3, i4, i5, i6, i7, i8, i9, matches1, matches2, matches3, matches4 },
#    mt = MT[ring];
#    non0Cons = NZSC[ring];
#    zeros = Position[ mt, x_Integer/; x === 0 ];
#    d = CC[ring] /@ Range[Rank[ring]];
#
#    Catch[
#      Do[
#        { i2, i1, i3 } = ind;
#        matches1 = Cases[ non0Cons, { i4_, i1, i6_ } ];
#        Do[
#          { i4, i6 } = ind1[[ { 1, 3 } ]];
#          matches2 = Cases[ non0Cons, { i5_, i4, i2 } ];
#          Do[
#            i5 = ind2[[ 1 ]];
#            If[
#              mt[[ i5, i6, i3 ]] =!= 0,
#              matches3 = Cases[ non0Cons, { i7_, i9_, i1 } ];
#              Do[
#                { i7, i9 } = ind3[[ { 1, 2 } ]];
#                matches4 =
#                  Cases[
#                    Range @ Rank @ ring, i0_ /;
#                    mt[[i4,i7,i0]] === 1 &&
#                    mt[[i6,d[[i9]],i0]] == 1 &&
#                    MemberQ[ crit1[ { i9, i0, i6, i7, i4, i1  }, d, mt ], 1 ]
#                  ];
#                Do[
#                  i0 = ind4;
#                  If[
#                    !MissingQ[
#                      FirstCase[
#                        Range @ Rank @ ring, i8_ /;
#                        mt[[i2,i7,i8]] =!= 0 && mt[[i8,i9,i3]] =!= 0 &&
#                        mt[[d[[i5]],i8,i0]] === 1 &&
#                        MemberQ[ crit1[ { i7, i2, i8, i4, i5, i0  }, d, mt ], 1 ] &&
#                        MemberQ[ crit1[ { i9, i8, i3, i0, i5, i6  }, d, mt ], 1 ] &&
#                        crit3[ { i4, i5, i6, i7, i8, i9 }, d, mt ] == 1
#                      ]
#                    ],
#                    Throw @ True
#                  ]
#                ,{ ind4, matches4 }]
#              ,{ ind3, matches3 }]
#            ]
#            ,{ ind2, matches2 }]
#          ,{ ind1, matches1 }]
#        ,{ ind, zeros }
#      ];
#      False
#    ]
#
#  ];

#TODO: to me: take a look at and see if warrants replacement: 
"""
    zsc_criterion(ring)

Return `true` if the Zero Spectrum Criterion rules out categorifiability.

This is a direct Julia port of the Mathematica `ZSCriterion` logic, but with
explicit Boolean helpers.
"""
function zsc_criterion(ring::FusionRing)::Bool
  mt = _mult_table(ring)
  r = size(mt, 1)

  d = _dual_indices_from_mt(mt)
  nonzero = _nonzero_structure_constants(mt)

  for i2 in 1:r, i1 in 1:r, i3 in 1:r
    mt[i2, i1, i3] == 1 || continue

    for ind1 in nonzero
      # Mathematica pattern: {i4_, i1, i6_}
      ind1[2] == i1 || continue

      i4 = ind1[1]
      i6 = ind1[3]

      for ind2 in nonzero
        # Mathematica pattern: {i5_, i4, i2}
        ind2[2] == i4 || continue
        ind2[3] == i2 || continue

        i5 = ind2[1]

        mt[i5, i6, i3] != 0 || continue
        _crit1_has_one((i1, i2, i3, i4, i5, i6), d, mt) || continue

        for ind3 in nonzero
          # Mathematica pattern: {i7_, i9_, i1}
          ind3[3] == i1 || continue

          i7 = ind3[1]
          i9 = ind3[2]

          for ind4 in nonzero
            # Mathematica pattern: {i2, i7, i8_}
            ind4[1] == i2 || continue
            ind4[2] == i7 || continue

            i8 = ind4[3]

            if mt[i8, i9, i3] != 0 &&
              _crit3_sum((i4, i5, i6, i7, i8, i9), d, mt) == 0 &&
              _crit2_has_one((i1, i2, i3, i7, i8, i9), d, mt)
              return true
            end
          end
        end
      end
    end
  end

  return false
end

#TODO: take a look at this too
"""
    osc_criterion(ring)

Return `true` if the One Spectrum Criterion rules out categorifiability.

This is a direct Julia port of the Mathematica `OSCriterion` logic.
"""
function osc_criterion(ring::FusionRing)::Bool
  mt = _mult_table(ring)
  r = size(mt, 1)

  d = _dual_indices_from_mt(mt)
  nonzero = _nonzero_structure_constants(mt)

  for i2 in 1:r, i1 in 1:r, i3 in 1:r
    mt[i2, i1, i3] == 0 || continue

    for ind1 in nonzero
      # Mathematica pattern: {i4_, i1, i6_}
      ind1[2] == i1 || continue

      i4 = ind1[1]
      i6 = ind1[3]

      for ind2 in nonzero
        # Mathematica pattern: {i5_, i4, i2}
        ind2[2] == i4 || continue
        ind2[3] == i2 || continue

        i5 = ind2[1]

        mt[i5, i6, i3] != 0 || continue

        for ind3 in nonzero
          # Mathematica pattern: {i7_, i9_, i1}
          ind3[3] == i1 || continue

          i7 = ind3[1]
          i9 = ind3[2]

          for i0 in 1:r
            mt[i4, i7, i0] == 1 || continue
            mt[i6, d[i9], i0] == 1 || continue
            _crit1_has_one((i9, i0, i6, i7, i4, i1), d, mt) || continue

            for i8 in 1:r
              if mt[i2, i7, i8] != 0 &&
                mt[i8, i9, i3] != 0 &&
                mt[d[i5], i8, i0] == 1 &&
                _crit1_has_one((i7, i2, i8, i4, i5, i0), d, mt) &&
                _crit1_has_one((i9, i8, i3, i0, i5, i6), d, mt) &&
                _crit3_sum((i4, i5, i6, i7, i8, i9), d, mt) == 1
                return true
              end
            end
          end
        end
      end
    end
  end

  return false
end
