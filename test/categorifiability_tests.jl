using Oscar

@testset "Categorifiability criteria" begin
  maybe_testset("basic", "1. public criteria return booleans") do
    z2 = zn_fusion_ring(2)

    check_true(csp_criterion(z2) isa Bool, "csp_criterion did not return a Bool")
    check_true(pdc_criterion(z2) isa Bool, "pdc_criterion did not return a Bool")
    check_true(dn_criterion(z2) isa Bool, "dn_criterion did not return a Bool")
    check_true(zsc_criterion(z2) isa Bool, "zsc_criterion did not return a Bool")
    return check_true(osc_criterion(z2) isa Bool, "osc_criterion did not return a Bool")
  end

  maybe_testset("basic", "crit results") do
  z2 = zn_fusion_ring(2)

  result = csp_result(z2)

  check_true(
    result isa CriterionResult,
    "csp_result did not return a CriterionResult",
  )

  check_equal(
    result.status,
    NoObstructionFound,
    "CSP returned an unexpected status for Z2",
  )

  check_equal(
    result.failed_input,
    nothing,
    "Unobstructed CSP result returned a failed input",
  )

  data = CategorifiabilityData()
  data[:csp] = result

  check_true(
    haskey(data, :csp),
    "CategorifiabilityData did not store the CSP result",
  )

  return check_equal(
    data[:csp],
    result,
    "CategorifiabilityData returned the wrong CSP result",
  )
end

  maybe_testset("intermediate", "2. formulas and cache behavior") do
    z2 = zn_fusion_ring(2)

    check_false(csp_criterion(z2), "CSP incorrectly obstructed the Z2 fusion ring")
    check_false(pdc_criterion(z2), "PDC incorrectly obstructed the Z2 fusion ring")
    check_false(dn_criterion(z2), "D-number criterion incorrectly obstructed Z2")
    check_false(zsc_criterion(z2), "ZSC incorrectly obstructed the Z2 fusion ring")
    check_false(osc_criterion(z2), "OSC incorrectly obstructed the Z2 fusion ring")

    half = QQBar(1) / QQBar(2)
    nond_polynomial = polynomial(ZZ, [2, 1, 1])
    nond = first(roots(QQBar, nond_polynomial))

    check_true(is_d_number(sqrt(QQBar(2))), "sqrt(2) was not recognized as a d-number")
    check_false(is_d_number(half), "is_d_number accepted a nonintegral rational")
    check_false(
      is_d_number(nond),
      "is_d_number used the minimal-polynomial coefficients in the wrong order",
    )

    minus_one_id = qqb_id(QQBar(-1))
    stale_characters = change_properties(
      z2, :characters => fill(minus_one_id, rank(z2), rank(z2))
    )
    check_true(
      csp_criterion(stale_characters),
      "CSP did not interpret character rows and simple-object columns correctly",
    )
    check_false(
      csp_criterion(stale_characters; force_compute = true),
      "CSP force_compute reused a stale character table",
    )

    stale_pdc = change_properties(
      z2, :formal_codegrees => qqb_id.([QQBar(2), QQBar(3)])
    )
    check_true(pdc_criterion(stale_pdc), "PDC did not use stored formal codegrees")
    check_false(
      pdc_criterion(stale_pdc; force_compute = true),
      "PDC force_compute reused stale formal codegrees",
    )

    stale_dn = change_properties(
      z2, :formal_codegrees => fill(qqb_id(nond), rank(z2))
    )
    check_true(dn_criterion(stale_dn), "D-number criterion ignored a non-d-number")
    return check_false(
      dn_criterion(stale_dn; force_compute = true),
      "D-number force_compute reused stale formal codegrees",
    )
  end

  maybe_testset("reference", "3. reference / Anyonica parity") do
    data = load_anyonica_data("categorifiability.json")

    for idx in oracle_case_indices("categorifiability.json", data; max_cases = 8)
      expected = data["Output"][idx]
      isnothing(expected) && continue

      ring = ring_from_anyonica_code(data["Input"][idx])
      actual = [zsc_criterion(ring), pdc_criterion(ring), dn_criterion(ring)]
      check_equal(
        actual,
        [expected[1], expected[3], expected[4]],
        "categorifiability criteria did not match Anyonica case $idx",
      )
    end
  end
end

#changed: added tests for pudc
#it should not matter where a value happens to appear in a vector!
@testset "pseudo-unitary Drinfeld-center criterion" begin
  ising = su2k_fusion_ring(2)

  result = FusionRings.pudc_criterion_explanation(
    ising;
    force_compute = true,
  )

  @test result isa Tuple
  @test length(result) == 3
  @test result[1] == false

  @test !FusionRings.pudc_criterion(
    ising;
    force_compute = true,
  )
end

#changed: testing independence from cached ordering
@testset "PUDC is independent of codegree ordering" begin
  ising = su2k_fusion_ring(2)

  codegrees = formal_codegrees(
    ising;
    force_compute = true,
  )

  normal_order = change_properties(
    ising,
    :formal_codegrees => qqb_id.(codegrees),
  )

  reversed_order = change_properties(
    ising,
    :formal_codegrees => qqb_id.(reverse(codegrees)),
  )

  @test FusionRings.pudc_criterion(normal_order) ==
        FusionRings.pudc_criterion(reversed_order)

  @test !FusionRings.pudc_criterion(normal_order)
  @test !FusionRings.pudc_criterion(reversed_order)
end