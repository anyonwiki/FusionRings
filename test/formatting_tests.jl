@testset "Formatting and printing" begin
  using Oscar

  # ============================================================
  # Integer formatting helpers
  # ============================================================

  @testset "integer formatting helpers" begin
    maybe_testset("basic", "1. basic construction") do
      b = FusionRings.bold_integer(12)
      s = FusionRings.subscript_integer(12)
      t = FusionRings.superscript_integer(12)

      check_true(b isa AbstractString, "FusionRings.bold_integer(12) did not return a string")
      check_true(s isa AbstractString, "FusionRings.subscript_integer(12) did not return a string")
      check_true(t isa AbstractString, "FusionRings.superscript_integer(12) did not return a string")

      check_false(isempty(b), "bold_integer(12) returned an empty string")
      check_false(isempty(s), "FusionRings.subscript_integer(12) returned an empty string")
      return check_false(isempty(t), "superscript_integer(12) returned an empty string")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      check_equal(
        FusionRings.transform_integer(0, FusionRings.bold_digits_dict),
        "𝟎",
        "FusionRings.transform_integer(0, FusionRings.bold_digits_dict) did not return bold zero",
      )
      check_equal(
        FusionRings.transform_integer(123, FusionRings.bold_digits_dict),
        "𝟏𝟐𝟑",
        "FusionRings.transform_integer(123, FusionRings.bold_digits_dict) did not return the expected bold digits",
      )

      check_equal(
        FusionRings.bold_integer(907),
        "𝟗𝟎𝟕",
        "FusionRings.bold_integer(907) did not return the expected bold representation",
      )
      check_equal(
        FusionRings.subscript_integer(314),
        "₃₁₄",
        "FusionRings.subscript_integer(314) did not return the expected subscript representation",
      )
      check_equal(
        FusionRings.superscript_integer(256),
        "²⁵⁶",
        "superscript_integer(256) did not return the expected superscript representation",
      )

      check_equal(
        FusionRings.bold_integer(1001),
        "𝟏𝟎𝟎𝟏",
        "bold_integer(1001) did not preserve repeated digits and zeros",
      )
      check_equal(
        FusionRings.subscript_integer(1001),
        "₁₀₀₁",
        "FusionRings.subscript_integer(1001) did not preserve repeated digits and zeros",
      )
      return check_equal(
        FusionRings.superscript_integer(1001),
        "¹⁰⁰¹",
        "FusionRings.superscript_integer(1001) did not preserve repeated digits and zeros",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # FusionRings.element_to_string
  # ============================================================

  @testset "FusionRings.element_to_string" begin
    maybe_testset("basic", "1. basic construction") do
      a = FusionRings.element_to_string(0, "x")
      b = FusionRings.element_to_string(1, "x")
      c = FusionRings.element_to_string(2, "x")

      check_true(
        a isa AbstractString, "FusionRings.element_to_string(0, \"x\") did not return a string"
      )
      check_true(
        b isa AbstractString, "FusionRings.element_to_string(1, \"x\") did not return a string"
      )
      return check_true(
        c isa AbstractString, "FusionRings.element_to_string(2, \"x\") did not return a string"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      check_equal(
        FusionRings.element_to_string(0, "x"),
        "",
        "FusionRings.element_to_string(0, \"x\") did not return the empty string",
      )
      check_equal(
        FusionRings.element_to_string(1, "x"), "x", "FusionRings.element_to_string(1, \"x\") did not return \"x\""
      )
      check_equal(
        FusionRings.element_to_string(2, "x"),
        "2 x",
        "FusionRings.element_to_string(2, \"x\") did not return \"2 x\"",
      )
      return check_equal(
        FusionRings.element_to_string(7, "abc"),
        "7 abc",
        "FusionRings.element_to_string(7, \"abc\") did not return \"7 abc\"",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # row_to_string
  # ============================================================

  @testset "row_to_string" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)
      s = FusionRings.row_to_string(z3, [0, 1, 0])

      return check_true(s isa AbstractString, "FusionRings.row_to_string did not return a string")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)

      check_equal(
        FusionRings.row_to_string(z3, [1, 0, 0]), "0", "FusionRings.row_to_string(Z3, [1,0,0]) did not return \"0\""
      )
      check_equal(
        FusionRings.row_to_string(z3, [0, 1, 0]), "1", "FusionRings.row_to_string(Z3, [0,1,0]) did not return \"1\""
      )
      check_equal(
        FusionRings.row_to_string(z3, [0, 0, 1]), "2", "FusionRings.row_to_string(Z3, [0,0,1]) did not return \"2\""
      )
      check_equal(
        FusionRings.row_to_string(z3, [1, 0, 1]),
        "0 ⊕ 2",
        "row_to_string(Z3, [1,0,1]) did not join summands with ⊕",
      )
      return check_equal(
        FusionRings.row_to_string(z3, [2, 0, 1]),
        "2 0 ⊕ 2",
        "FusionRings.row_to_string(Z3, [2,0,1]) did not format multiplicities correctly",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # FusionRings.product_string
  # ============================================================

  @testset "FusionRings.product_string" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)
      s = FusionRings.product_string(z2, 1, 1)

      check_true(s isa AbstractString, "FusionRings.product_string(Z2, 1, 1) did not return a string")
      return check_false(isempty(s), "FusionRings.product_string(Z2, 1, 1) returned an empty string")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)
      su2_2 = su2k_fusion_ring(2)

      check_equal(
        FusionRings.product_string(z2, 1, 1),
        "0 × 0 = 0",
        "FusionRings.product_string(Z2, 1, 1) was not \"0 × 0 = 0\"",
      )
      check_equal(
        FusionRings.product_string(z2, 2, 2),
        "1 × 1 = 0",
        "FusionRings.product_string(Z2, 2, 2) was not \"1 × 1 = 0\"",
      )
      check_equal(
        FusionRings.product_string(z3, 2, 2),
        "1 × 1 = 2",
        "FusionRings.product_string(Z3, 2, 2) was not \"1 × 1 = 2\"",
      )
      check_equal(
        FusionRings.product_string(z3, 2, 3),
        "1 × 2 = 0",
        "FusionRings.product_string(Z3, 2, 3) was not \"1 × 2 = 0\"",
      )
      return check_equal(
        FusionRings.product_string(su2_2, 2, 2),
        "1 × 1 = 0 ⊕ 2",
        "FusionRings.product_string(SU(2)_2, 2, 2) did not show both summands correctly",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # print_multiplication_table / pmt
  # ============================================================

  @testset "print_multiplication_table / pmt" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)

      tab = print_multiplication_table(z2)
      tab_alias = pmt(z2)

      check_true(
        tab isa AbstractMatrix,
        "print_multiplication_table(Z2) did not return a matrix-like object",
      )
      check_equal(
        size(tab),
        (rank(z2), rank(z2)),
        "print_multiplication_table(Z2) did not return a rank×rank table",
      )
      return check_equal(
        size(tab_alias), (rank(z2), rank(z2)), "pmt(Z2) did not return a rank×rank table"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)

      tab2 = print_multiplication_table(z2)
      check_equal(tab2[1, 1], "0", "print_multiplication_table(Z2)[1,1] was not \"0\"")
      check_equal(tab2[2, 2], "0", "print_multiplication_table(Z2)[2,2] was not \"0\"")
      check_equal(tab2[1, 2], "1", "print_multiplication_table(Z2)[1,2] was not \"1\"")

      tab3 = print_multiplication_table(z3)
      check_equal(tab3[2, 2], "2", "print_multiplication_table(Z3)[2,2] was not \"2\"")
      check_equal(tab3[2, 3], "0", "print_multiplication_table(Z3)[2,3] was not \"0\"")

      return check_equal(
        pmt(z2), tab2, "pmt(z2) did not agree with print_multiplication_table(z2)"
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # Base.show(::FusionRing)
  # ============================================================

  @testset "Base.show(::FusionRing)" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)
      shown = sprint(show, z2)

      check_true(
        shown isa AbstractString, "show(zn_fusion_ring(2)) did not produce a string"
      )
      return check_false(isempty(shown), "show(zn_fusion_ring(2)) produced an empty string")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      named = zn_fusion_ring(2)
      shown_named = sprint(show, named)

      check_true(
        startswith(shown_named, "FR("), "show(zn_fusion_ring(2)) did not begin with \"FR(\""
      )
      check_true(
        occursin("ℤ₂", shown_named),
        "show(zn_fusion_ring(2)) did not contain the ring name",
      )
      check_true(
        endswith(shown_named, ")"), "show(zn_fusion_ring(2)) did not end with \")\""
      )

      mt = zeros(Int, 2, 2, 2)
      mt[1, 1, 1] = 1
      mt[1, 2, 2] = 1
      mt[2, 1, 2] = 1
      mt[2, 2, 1] = 1

      unnamed = fusion_ring(mt; labels = ["a", "b"])
      shown_unnamed = sprint(show, unnamed)

      check_true(
        startswith(shown_unnamed, "FR("),
        "show on an unnamed fusion ring did not begin with \"FR(\"",
      )
      return check_true(
        endswith(shown_unnamed, ")"),
        "show on an unnamed fusion ring did not end with \")\"",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # string cleanup helpers
  # ============================================================

  @testset "string cleanup helpers" begin
    maybe_testset("basic", "1. basic construction") do
      a = FusionRings.fix_fractions("3//4")
      b = FusionRings.fix_mult("2*x")
      c = FusionRings.fix_spaces("a b")
      d = FusionRings.fix_powers("x^12")
      e = FusionRings.fix_cyclo("zeta(5)")
      f = FusionRings.fix_poly_string(" 3//4 * x^12 ")

      check_true(a isa AbstractString, "FusionRings.fix_fractions did not return a string")
      check_true(b isa AbstractString, "FusionRings.FusionRings.fix_mult did not return a string")
      check_true(c isa AbstractString, "FusionRings.fix_spaces did not return a string")
      check_true(d isa AbstractString, "FusionRings.fix_powers did not return a string")
      check_true(e isa AbstractString, "FusionRings.fix_cyclo did not return a string")
      return check_true(f isa AbstractString, "fix_poly_string did not return a string")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      check_equal(
        FusionRings.fix_fractions("3//4"),
        "\\frac{3}{4}",
        "FusionRings.fix_fractions(\"3//4\") did not convert a simple rational",
      )
      check_equal(
        FusionRings.fix_fractions("1//2 + 3//4"),
        "\\frac{1}{2} + \\frac{3}{4}",
        "FusionRings.fix_fractions did not convert multiple rationals in one string",
      )

      check_equal(FusionRings.fix_mult("2*x"), "2x", "FusionRings.fix_mult(\"2*x\") did not remove *")
      check_equal(FusionRings.fix_mult("a*b*c"), "abc", "FusionRings.fix_mult(\"a*b*c\") did not remove all *")

      check_equal(
        FusionRings.fix_spaces("a b  c"), "abc", "FusionRings.fix_spaces(\"a b  c\") did not remove spaces"
      )
      check_equal(
        FusionRings.fix_spaces("   x   "),
        "x",
        "FusionRings.fix_spaces(\"   x   \") did not remove leading/trailing spaces",
      )

      check_equal(
        FusionRings.fix_powers("x^12"),
        "x^{12}",
        "FusionRings.fix_powers(\"x^12\") did not brace a multi-digit exponent",
      )
      check_equal(
        FusionRings.fix_powers("x^2"),
        "x^2",
        "FusionRings.fix_powers(\"x^2\") should not change a one-digit exponent",
      )
      check_equal(
        FusionRings.fix_powers("x^12 + y^34"),
        "x^{12} + y^{34}",
        "FusionRings.fix_powers did not handle multiple multi-digit exponents",
      )

      check_equal(
        FusionRings.fix_cyclo("zeta(5)"),
        "\\zeta_{5}",
        "FusionRings.fix_cyclo(\"zeta(5)\") did not convert to zeta subscript notation",
      )
      check_equal(
        FusionRings.fix_cyclo("zeta(3) + zeta(7)"),
        "\\zeta_{3} + \\zeta_{7}",
        "FusionRings.fix_cyclo did not handle multiple zeta terms",
      )

      return check_equal(
        FusionRings.fix_poly_string(" 3//4 * x^12 "),
        "\\frac{3}{4}x^{12}",
        "FusionRings.fix_poly_string did not compose cleanup helpers correctly",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # factor_squares / is_geometric_array / is_power_sum
  # ============================================================

  @testset "small algebraic formatting helpers" begin
    maybe_testset("basic", "1. basic construction") do
      fs = factor_squares(12)
      g1 = FusionRings.is_geometric_array([1, 2, 4, 8])

      check_true(fs isa Tuple, "factor_squares(12) did not return a tuple")
      return check_true(g1 isa Bool, "FusionRings.is_geometric_array([1,2,4,8]) did not return a Bool")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      check_equal(factor_squares(1), (1, 1), "factor_squares(1) was not (1,1)")
      check_equal(factor_squares(7), (1, 7), "factor_squares(7) was not (1,7)")
      check_equal(factor_squares(12), (2, 3), "factor_squares(12) was not (2,3)")
      check_equal(factor_squares(18), (3, 2), "factor_squares(18) was not (3,2)")
      check_equal(factor_squares(16), (4, 1), "factor_squares(16) was not (4,1)")

      check_true(FusionRings.is_geometric_array(Int[]), "FusionRings.is_geometric_array(Int[]) should be true")
      check_true(FusionRings.is_geometric_array([1]), "FusionRings.is_geometric_array([1]) should be true")
      check_false(FusionRings.is_geometric_array([2]), "FusionRings.is_geometric_array([2]) should be false")
      check_true(
        FusionRings.is_geometric_array([1, 2, 4, 8]), "FusionRings.is_geometric_array([1,2,4,8]) should be true"
      )
      check_true(
        FusionRings.is_geometric_array([1, -1, 1, -1]), "FusionRings.is_geometric_array([1,-1,1,-1]) should be true"
      )
      check_false(
        FusionRings.is_geometric_array([1, 2, 5, 10]), "FusionRings.is_geometric_array([1,2,5,10]) should be false"
      )

      # is_power_sum: x^2 + x + 1 has coefficients [1,1,1]
      R, x = Oscar.polynomial_ring(QQ, "x")
      check_true(FusionRings.is_power_sum(x^2 + x + 1), "is_power_sum(x^2 + x + 1) should be true")
      return check_false(
        FusionRings.is_power_sum(x^2 + 2*x + 1), "is_power_sum(x^2 + 2x + 1) should be false"
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end

  # ============================================================
  # integration: default labels
  # ============================================================

  @testset "default labels use bold_integer" begin
    maybe_testset("basic", "1. basic construction") do
      mt = zeros(Int, 2, 2, 2)
      mt[1, 1, 1] = 1
      mt[1, 2, 2] = 1
      mt[2, 1, 2] = 1
      mt[2, 2, 1] = 1

      r = fusion_ring(mt)

      check_true(r isa FusionRing, "fusion_ring(mt) did not return a FusionRing")
      return check_equal(
        length(labels(r)),
        rank(r),
        "fusion_ring(mt) did not produce the correct number of labels",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      mt = zeros(Int, 2, 2, 2)
      mt[1, 1, 1] = 1
      mt[1, 2, 2] = 1
      mt[2, 1, 2] = 1
      mt[2, 2, 1] = 1

      r = fusion_ring(mt)

      check_equal(
        labels(r)[1],
        FusionRings.bold_integer(1),
        "fusion_ring(mt) did not use bold_integer(1) as the first default label",
      )
      return check_equal(
        labels(r)[2],
        FusionRings.bold_integer(2),
        "fusion_ring(mt) did not use bold_integer(2) as the second default label",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # TODO
      @test true
    end
  end
end
