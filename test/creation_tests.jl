@testset "Constructors" begin
  """
      make_z2_mt()

  Multiplication table for Z2 with basis:
  1 = vacuum, 2 = nontrivial element, and 2⊗2 = 1.
  """
  function make_z2_mt()
    mt = zeros(Int, 2, 2, 2)
    mt[1, 1, 1] = 1
    mt[1, 2, 2] = 1
    mt[2, 1, 2] = 1
    mt[2, 2, 1] = 1
    return mt
  end

  """
      make_z3_mt()

  Multiplication table for Z3 with basis:
  1↔0, 2↔1, 3↔2 mod 3.
  """
  function make_z3_mt()
    mt = zeros(Int, 3, 3, 3)
    for i in 0:2, j in 0:2
      k = mod(i + j, 3)
      mt[i + 1, j + 1, k + 1] = 1
    end
    return mt
  end

  # ============================================================
  # fusion_ring core constructor
  # ============================================================

  @testset "fusion_ring core constructor" begin
    maybe_testset("basic", "1. basic construction") do
      mt = make_z2_mt()
      r = fusion_ring(mt; labels = ["0", "1"], names = ["Z2"])

      check_true(r isa FusionRing, "fusion_ring did not return a FusionRing")
      check_true(rank(r) > 0, "fusion_ring did not construct a positive-rank ring")
      check_equal(
        size(multiplication_table(r)),
        (rank(r), rank(r), rank(r)),
        "fusion_ring did not produce a cubic multiplication table",
      )
      check_equal(
        length(labels(r)),
        rank(r),
        "fusion_ring did not produce the correct number of labels",
      )

      r_default = fusion_ring(mt)
      check_true(
        r_default isa FusionRing,
        "fusion_ring without explicit labels did not return a FusionRing",
      )
      check_equal(
        size(multiplication_table(r_default)),
        (rank(r_default), rank(r_default), rank(r_default)),
        "fusion_ring without explicit labels did not produce a cubic multiplication table",
      )
      check_equal(
        length(labels(r_default)),
        rank(r_default),
        "fusion_ring without explicit labels did not produce the correct number of labels",
      )

      mt_bad = zeros(Int, 2, 2, 2)
      r_skip = fusion_ring(mt_bad; skip_check = true)
      check_true(
        r_skip isa FusionRing,
        "fusion_ring(...; skip_check=true) did not return a FusionRing",
      )
      return check_equal(
        size(multiplication_table(r_skip)),
        (rank(r_skip), rank(r_skip), rank(r_skip)),
        "fusion_ring(...; skip_check=true) did not preserve a cubic tensor shape",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      mt = make_z2_mt()
      r = fusion_ring(mt; labels = ["0", "1"], names = ["Z2"])

      check_equal(
        rank(r),
        2,
        "fusion_ring did not construct a rank-2 ring from a valid Z2 multiplication table",
      )
      check_equal(
        labels(r), ["0", "1"], "fusion_ring did not preserve labels for a valid Z2 table"
      )
      check_equal(
        names(r), ["Z2"], "fusion_ring did not preserve names for a valid Z2 table"
      )
      check_true(
        is_commutative(r),
        "fusion_ring constructed a valid Z2 table but the result was not detected as commutative",
      )
      check_true(
        is_group_ring(r),
        "fusion_ring constructed a valid Z2 table but the result was not detected as a group ring",
      )


      # Fix: added Keep  same total number of unit coefficients while assigning
# the wrong duals.  guards against checking only
# sum(mt[:, :, 1]) == rank.
      mt_bad_dual_distribution = make_z3_mt()
      mt_bad_dual_distribution[2, 3, 1] = 0
      mt_bad_dual_distribution[2, 2, 1] = 1

      check_false(
      FusionRings.check_inverse(mt_bad_dual_distribution),
      "check_inverse accepted an invalid distribution of dual objects",)


      check_equal(
        FusionRings.fusion_product(r, 1, 1),
        Dict(1 => 1),
        "fusion_ring(valid Z2): vacuum × vacuum was not vacuum",
      )
      return check_equal(
        FusionRings.fusion_product(r, 2, 2),
        Dict(1 => 1),
        "fusion_ring(valid Z2): nontrivial simple squared was not vacuum",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("zn_tables.json")

      for idx in oracle_case_indices("zn_tables.json", data)
        mt = json_int_3tensor(data["Output"][idx])
        r = fusion_ring(mt; skip_check = true)

        check_equal_tensor(
          multiplication_table(r),
          mt,
          "fusion_ring could not reconstruct Anyonica zn_tables.json table case $idx",
        )
      end
    end

    if want_level("basic") || want_level("intermediate")
      @testset "validation failures" begin
        mt = make_z2_mt()
        mt[2, 2, 1] = -1
        check_throws(
          () -> fusion_ring(mt; labels = ["0", "1"]),
          "fusion_ring accepted a multiplication table with a negative structure constant",
        )

        mtf = Array{Float64}(undef, 2, 2, 2)
        fill!(mtf, 0.0)
        mtf[1, 1, 1] = 1.0
        mtf[1, 2, 2] = 1.0
        mtf[2, 1, 2] = 1.0
        mtf[2, 2, 1] = 1.0
        check_throws(
          () -> fusion_ring(mtf; labels = ["0", "1"]),
          "fusion_ring accepted a multiplication table with non-integer structure constants",
        )

        mt_bad_dims = zeros(Int, 2, 2, 3)
        check_throws(
          () -> fusion_ring(mt_bad_dims),
          "fusion_ring accepted a multiplication table whose tensor dimensions were not all equal",
        )

        mt_bad_unit = make_z2_mt()
        mt_bad_unit[1, 2, 2] = 0
        check_throws(
          () -> fusion_ring(mt_bad_unit; labels = ["0", "1"]),
          "fusion_ring accepted a multiplication table whose first basis element was not a unit",
        )

        mt_bad_inverse = make_z2_mt()
        mt_bad_inverse[2, 1, 1] = 1
        check_throws(
          () -> fusion_ring(mt_bad_inverse; labels = ["0", "1"]),
          "fusion_ring accepted a multiplication table violating the unique inverse condition",
        )

        mt_bad_assoc = zeros(Int, 2, 2, 2)
        mt_bad_assoc[1, 1, 1] = 1
        mt_bad_assoc[1, 2, 2] = 1
        mt_bad_assoc[2, 1, 2] = 1
        mt_bad_assoc[2, 2, 2] = 1
        check_throws(
          () -> fusion_ring(mt_bad_assoc; labels = ["0", "1"]),
          "fusion_ring accepted a non-fusion multiplication table that should fail validation",
        )

        mt_good = make_z2_mt()
        check_throws(
          () -> fusion_ring(mt_good; labels = ["0"]),
          "fusion_ring accepted labels whose length did not equal the rank",
        )
      end
    end
  end

  # ============================================================
  # zn_fusion_ring
  # ============================================================

  @testset "zn_fusion_ring" begin
    maybe_testset("basic", "1. basic construction") do
      for n in 1:4
        z = zn_fusion_ring(n)
        check_true(z isa FusionRing, "zn_fusion_ring($n) did not return a FusionRing")
        check_true(rank(z) > 0, "zn_fusion_ring($n) did not produce a positive-rank ring")
        check_equal(
          size(multiplication_table(z)),
          (rank(z), rank(z), rank(z)),
          "zn_fusion_ring($n) did not produce a cubic multiplication table",
        )
        check_equal(
          length(labels(z)),
          rank(z),
          "zn_fusion_ring($n) did not produce the correct number of labels",
        )
      end
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z1 = zn_fusion_ring(1)
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)

      check_equal(rank(z1), 1, "rank(zn_fusion_ring(1)) was not 1")
      check_equal(rank(z2), 2, "rank(zn_fusion_ring(2)) was not 2")
      check_equal(rank(z3), 3, "rank(zn_fusion_ring(3)) was not 3")
      check_equal(rank(z4), 4, "rank(zn_fusion_ring(4)) was not 4")

      check_equal(labels(z1), ["0"], "labels(zn_fusion_ring(1)) were incorrect")
      check_equal(labels(z2), ["0", "1"], "labels(zn_fusion_ring(2)) were incorrect")
      check_equal(labels(z3), ["0", "1", "2"], "labels(zn_fusion_ring(3)) were incorrect")
      check_equal(
        labels(z4), ["0", "1", "2", "3"], "labels(zn_fusion_ring(4)) were incorrect"
      )

      check_true(is_commutative(z1), "zn_fusion_ring(1) was not commutative")
      check_true(is_commutative(z2), "zn_fusion_ring(2) was not commutative")
      check_true(is_commutative(z3), "zn_fusion_ring(3) was not commutative")
      check_true(is_commutative(z4), "zn_fusion_ring(4) was not commutative")

      check_true(is_group_ring(z1), "zn_fusion_ring(1) was not detected as a group ring")
      check_true(is_group_ring(z2), "zn_fusion_ring(2) was not detected as a group ring")
      check_true(is_group_ring(z3), "zn_fusion_ring(3) was not detected as a group ring")
      check_true(is_group_ring(z4), "zn_fusion_ring(4) was not detected as a group ring")

      mt1 = multiplication_table(z3)
      check_equal(mt1[1, 1, 1], 1, "vacuum × vacuum was not vacuum in zn_fusion_ring(3)")
      check_equal(mt1[2, 2, 3], 1, "index 2 × index 2 in zn_fusion_ring(3) was not index 3")
      check_equal(mt1[2, 3, 1], 1, "index 2 × index 3 in zn_fusion_ring(3) was not vacuum")

      mt2 = multiplication_table(z4)
      return check_equal(
        mt2[3, 3, 1], 1, "index 3 × index 3 in zn_fusion_ring(4) was not vacuum"
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("zn_tables.json")

      for idx in oracle_case_indices("zn_tables.json", data)
        n = json_int(data["Input"][idx])

        check_mt_equal(
          zn_fusion_ring(n),
          data["Output"][idx],
          "zn_fusion_ring($n) did not match Anyonica zn_tables.json case $idx",
        )
      end
    end
  end

  # ============================================================
  # su2k_fusion_ring
  # ============================================================

  @testset "su2k_fusion_ring" begin
    maybe_testset("basic", "1. basic construction") do
      for k in 1:3
        r = su2k_fusion_ring(k)
        check_true(r isa FusionRing, "su2k_fusion_ring($k) did not return a FusionRing")
        check_true(rank(r) > 0, "su2k_fusion_ring($k) did not produce a positive-rank ring")
        check_equal(
          size(multiplication_table(r)),
          (rank(r), rank(r), rank(r)),
          "su2k_fusion_ring($k) did not produce a cubic multiplication table",
        )
        check_equal(
          length(labels(r)),
          rank(r),
          "su2k_fusion_ring($k) did not produce the correct number of labels",
        )
      end
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      r1 = su2k_fusion_ring(1)
      r2 = su2k_fusion_ring(2)
      r3 = su2k_fusion_ring(3)

      check_equal(rank(r1), 2, "rank(su2k_fusion_ring(1)) was not 2")
      check_equal(rank(r2), 3, "rank(su2k_fusion_ring(2)) was not 3")
      check_equal(rank(r3), 4, "rank(su2k_fusion_ring(3)) was not 4")

      check_equal(labels(r1), ["0", "1"], "labels(su2k_fusion_ring(1)) were incorrect")
      check_equal(labels(r2), ["0", "1", "2"], "labels(su2k_fusion_ring(2)) were incorrect")
      check_equal(
        labels(r3), ["0", "1", "2", "3"], "labels(su2k_fusion_ring(3)) were incorrect"
      )

      check_true(is_commutative(r1), "su2k_fusion_ring(1) was not commutative")
      check_true(is_commutative(r2), "su2k_fusion_ring(2) was not commutative")
      check_true(is_commutative(r3), "su2k_fusion_ring(3) was not commutative")

      mt1 = multiplication_table(r1)
      check_equal(
        mt1[2, 2, 1], 1, "nontrivial simple squared in su2k_fusion_ring(1) was not vacuum"
      )

      mt2 = multiplication_table(r2)
      return check_equal(
        mt2[2, 2, :],
        [1, 0, 1],
        "index 2 × index 2 in su2k_fusion_ring(2) was not vacuum + top object",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("su2k_tables.json")

      for idx in oracle_case_indices("su2k_tables.json", data)
        k = json_int(data["Input"][idx])

        check_mt_equal(
          su2k_fusion_ring(k),
          data["Output"][idx],
          "su2k_fusion_ring($k) did not match Anyonica su2k_tables.json case $idx",
        )
      end
    end
  end

  # ============================================================
  # psu2k_fusion_ring
  # ============================================================

  @testset "psu2k_fusion_ring" begin
    maybe_testset("basic", "1. basic construction") do
      for k in (2, 4, 6)
        r = psu2k_fusion_ring(k)
        check_true(r isa FusionRing, "psu2k_fusion_ring($k) did not return a FusionRing")
        check_true(
          rank(r) > 0, "psu2k_fusion_ring($k) did not produce a positive-rank ring"
        )
        check_equal(
          size(multiplication_table(r)),
          (rank(r), rank(r), rank(r)),
          "psu2k_fusion_ring($k) did not produce a cubic multiplication table",
        )
      end
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      r2 = psu2k_fusion_ring(2)
      r4 = psu2k_fusion_ring(4)
      r6 = psu2k_fusion_ring(6)

      check_equal(rank(r2), 2, "rank(psu2k_fusion_ring(2)) was not 2")
      check_equal(rank(r4), 3, "rank(psu2k_fusion_ring(4)) was not 3")
      check_equal(rank(r6), 4, "rank(psu2k_fusion_ring(6)) was not 4")

      check_true(is_commutative(r2), "psu2k_fusion_ring(2) was not commutative")
      check_true(is_commutative(r4), "psu2k_fusion_ring(4) was not commutative")
      check_true(is_commutative(r6), "psu2k_fusion_ring(6) was not commutative")

      mt1 = multiplication_table(r2)
      return check_equal(
        mt1[2, 2, 1], 1, "nontrivial simple squared in psu2k_fusion_ring(2) was not vacuum"
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("psu2k_tables.json")

      for idx in oracle_case_indices("psu2k_tables.json", data)
        k = json_int(data["Input"][idx])

        check_mt_equal(
          psu2k_fusion_ring(k),
          data["Output"][idx],
          "psu2k_fusion_ring($k) did not match Anyonica psu2k_tables.json case $idx",
        )
      end
    end
  end

  # ============================================================
  # son2_fusion_ring / metaplectic_fusion_ring
  # ============================================================

  @testset "son2_fusion_ring / metaplectic_fusion_ring" begin
    maybe_testset("basic", "1. basic construction") do
      for m in (5, 6, 8)
        r = son2_fusion_ring(m)
        check_true(r isa FusionRing, "son2_fusion_ring($m) did not return a FusionRing")
        check_true(rank(r) > 0, "son2_fusion_ring($m) did not produce a positive-rank ring")
        check_equal(
          size(multiplication_table(r)),
          (rank(r), rank(r), rank(r)),
          "son2_fusion_ring($m) did not produce a cubic multiplication table",
        )
        check_equal(
          length(labels(r)),
          rank(r),
          "son2_fusion_ring($m) did not produce the correct number of labels",
        )
      end

      meta = metaplectic_fusion_ring(5)
      check_true(
        meta isa FusionRing, "metaplectic_fusion_ring(5) did not return a FusionRing"
      )
      return check_true(
        rank(meta) > 0, "metaplectic_fusion_ring(5) did not produce a positive-rank ring"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      odd = son2_fusion_ring(5)
      even2 = son2_fusion_ring(6)
      even4 = son2_fusion_ring(8)
      meta = metaplectic_fusion_ring(5)

      check_true(is_commutative(odd), "son2_fusion_ring(5) was not commutative")
      check_true(is_commutative(even2), "son2_fusion_ring(6) was not commutative")
      check_true(is_commutative(even4), "son2_fusion_ring(8) was not commutative")

      check_equal(
        rank(meta),
        rank(odd),
        "metaplectic_fusion_ring(5) did not have the same rank as son2_fusion_ring(5)",
      )
      return check_equal_tensor(
        multiplication_table(meta),
        multiplication_table(odd),
        "metaplectic_fusion_ring(5) did not match son2_fusion_ring(5)",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("son2_tables.json")

      for idx in oracle_case_indices("son2_tables.json", data)
        N = json_int(data["Input"][idx])
        expected = data["Output"][idx]

        if N < 4
          # Anyonica exports SO(N)_2 for N = 1,2,3, but the
          # current Julia constructor intentionally requires N ≥ 4.
          @test_broken multiplication_table(son2_fusion_ring(N)) ==
            json_int_3tensor(expected)
        else
          check_mt_equal(
            son2_fusion_ring(N),
            expected,
            "son2_fusion_ring($N) did not match Anyonica son2_tables.json case $idx",
          )
          check_mt_equal(
            metaplectic_fusion_ring(N),
            expected,
            "metaplectic_fusion_ring($N) did not match Anyonica son2_tables.json case $idx",
          )
        end
      end
    end

    if want_level("basic") || want_level("intermediate")
      @testset "invalid input" begin
        check_throws(
          () -> son2_fusion_ring(3),
          "son2_fusion_ring accepted m=3 even though it should require m ≥ 4",
        )
      end
    end
  end

  # ============================================================
  # group_fusion_ring
  # ============================================================

  @testset "group_fusion_ring" begin
    z2_tab = [
      1 2
      2 1
    ]

    z3_tab = [
      1 2 3
      2 3 1
      3 1 2
    ]

    bad_not_group = [
      1 1
      2 2
    ]

    maybe_testset("basic", "1. basic construction") do
      r2 = group_fusion_ring(z2_tab)
      r3 = group_fusion_ring(z3_tab)

      check_true(
        r2 isa FusionRing, "group_fusion_ring(Z2 table) did not return a FusionRing"
      )
      check_true(
        r3 isa FusionRing, "group_fusion_ring(Z3 table) did not return a FusionRing"
      )

      check_equal(
        size(multiplication_table(r2)),
        (rank(r2), rank(r2), rank(r2)),
        "group_fusion_ring(Z2 table) did not produce a cubic multiplication table",
      )
      check_equal(
        size(multiplication_table(r3)),
        (rank(r3), rank(r3), rank(r3)),
        "group_fusion_ring(Z3 table) did not produce a cubic multiplication table",
      )
      check_equal(
        length(labels(r2)),
        rank(r2),
        "group_fusion_ring(Z2 table) did not produce the correct number of labels",
      )
      return check_equal(
        length(labels(r3)),
        rank(r3),
        "group_fusion_ring(Z3 table) did not produce the correct number of labels",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      r2 = group_fusion_ring(z2_tab)
      r3 = group_fusion_ring(z3_tab)

      check_equal(
        rank(r2), 2, "group_fusion_ring on the Z2 Cayley table did not produce rank 2"
      )
      check_equal(
        rank(r3), 3, "group_fusion_ring on the Z3 Cayley table did not produce rank 3"
      )

      check_equal(
        size(multiplication_table(r2)),
        (2, 2, 2),
        "group_fusion_ring on the Z2 Cayley table had wrong multiplication-table size",
      )
      check_equal(
        size(multiplication_table(r3)),
        (3, 3, 3),
        "group_fusion_ring on the Z3 Cayley table had wrong multiplication-table size",
      )

      check_true(
        is_group_ring(r2),
        "group_fusion_ring on the Z2 Cayley table was not detected as a group ring",
      )
      check_true(
        is_group_ring(r3),
        "group_fusion_ring on the Z3 Cayley table was not detected as a group ring",
      )
      check_true(
        is_commutative(r2), "group_fusion_ring on the Z2 Cayley table was not commutative"
      )
      check_true(
        is_commutative(r3), "group_fusion_ring on the Z3 Cayley table was not commutative"
      )

      check_equal_tensor(
        multiplication_table(r2),
        multiplication_table(zn_fusion_ring(2)),
        "group_fusion_ring(Z2 table) did not match zn_fusion_ring(2)",
      )
      return check_equal_tensor(
        multiplication_table(r3),
        multiplication_table(zn_fusion_ring(3)),
        "group_fusion_ring(Z3 table) did not match zn_fusion_ring(3)",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("grouptables.json")

      for idx in oracle_case_indices("grouptables.json", data)
        tab = json_int_matrix(data["Input"][idx])

        check_mt_equal(
          group_fusion_ring(tab),
          data["Output"][idx],
          "group_fusion_ring did not match Anyonica grouptables.json case $idx",
        )
      end
    end

    if want_level("basic") || want_level("intermediate")
      @testset "invalid input" begin
        check_throws(
          () -> group_fusion_ring(bad_not_group),
          "group_fusion_ring accepted a table that was not a valid Cayley table",
        )
      end
    end
  end

  # ============================================================
  # internal group-table helpers
  # ============================================================

  @testset "_is_group_table / _inverse_vector" begin
    z2_tab = [
      1 2
      2 1
    ]

    z3_tab = [
      1 2 3
      2 3 1
      3 1 2
    ]

    bad_tab = [
      1 2
      1 2
    ]

    maybe_testset("basic", "1. basic helper behaviour") do
      check_true(
        FusionRings._is_group_table(z2_tab),
        "_is_group_table did not recognize the Z2 Cayley table",
      )
      check_true(
        FusionRings._is_group_table(z3_tab),
        "_is_group_table did not recognize the Z3 Cayley table",
      )
      return check_false(
        FusionRings._is_group_table(bad_tab),
        "_is_group_table incorrectly accepted an invalid table",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      check_equal(
        FusionRings._inverse_vector(z2_tab),
        [1, 2],
        "_inverse_vector(Z2 table) was not [1,2]",
      )
      return check_equal(
        FusionRings._inverse_vector(z3_tab),
        [1, 3, 2],
        "_inverse_vector(Z3 table) was not [1,3,2]",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("grouptables.json")

      for idx in oracle_case_indices("grouptables.json", data)
        tab = json_int_matrix(data["Input"][idx])

        check_true(
          FusionRings._is_group_table(tab),
          "_is_group_table rejected Anyonica grouptables.json input case $idx",
        )

        inv = FusionRings._inverse_vector(tab)
        check_true(
          all(i -> tab[i, inv[i]] == 1 && tab[inv[i], i] == 1, 1:size(tab, 1)),
          "_inverse_vector did not return two-sided inverses for Anyonica grouptables.json case $idx",
        )
      end
    end
  end

  # ============================================================
  # TY_fusion_ring
  # ============================================================

  @testset "TY_fusion_ring" begin
    z2_tab = [
      1 2
      2 1
    ]

    z3_tab = [
      1 2 3
      2 3 1
      3 1 2
    ]

    bad_tab = [
      1 1
      2 2
    ]

    maybe_testset("basic", "1. basic construction") do
      r2 = TY_fusion_ring(z2_tab)
      r3 = TY_fusion_ring(z3_tab)

      check_true(r2 isa FusionRing, "TY_fusion_ring(Z2 table) did not return a FusionRing")
      check_true(r3 isa FusionRing, "TY_fusion_ring(Z3 table) did not return a FusionRing")

      check_equal(
        size(multiplication_table(r2)),
        (rank(r2), rank(r2), rank(r2)),
        "TY_fusion_ring(Z2) did not produce a cubic multiplication table",
      )
      check_equal(
        size(multiplication_table(r3)),
        (rank(r3), rank(r3), rank(r3)),
        "TY_fusion_ring(Z3) did not produce a cubic multiplication table",
      )
      check_equal(
        length(labels(r2)),
        rank(r2),
        "TY_fusion_ring(Z2) did not produce the correct number of labels",
      )
      return check_equal(
        length(labels(r3)),
        rank(r3),
        "TY_fusion_ring(Z3) did not produce the correct number of labels",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      r2 = TY_fusion_ring(z2_tab)
      r3 = TY_fusion_ring(z3_tab)

      check_equal(
        rank(r2), 3, "TY_fusion_ring on the Z2 Cayley table did not produce rank 3"
      )
      check_equal(
        rank(r3), 4, "TY_fusion_ring on the Z3 Cayley table did not produce rank 4"
      )

      check_equal(labels(r2), ["1", "2", "m"], "TY_fusion_ring(Z2) labels were incorrect")
      check_equal(
        labels(r3), ["1", "2", "3", "m"], "TY_fusion_ring(Z3) labels were incorrect"
      )

      check_equal(
        size(multiplication_table(r2)),
        (3, 3, 3),
        "TY_fusion_ring(Z2) multiplication table had wrong size",
      )
      check_equal(
        size(multiplication_table(r3)),
        (4, 4, 4),
        "TY_fusion_ring(Z3) multiplication table had wrong size",
      )

      m = rank(r2)
      mt1 = multiplication_table(r2)

      check_equal(mt1[1, m, m], 1, "vacuum × m was not m in TY_fusion_ring(Z2)")
      check_equal(
        mt1[2, m, m], 1, "nontrivial group element × m was not m in TY_fusion_ring(Z2)"
      )

      return check_equal(
        mt1[m, m, :],
        [1, 1, 0],
        "m × m in TY_fusion_ring(Z2) was not the sum of all group elements",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("TY_tables.json")

      for idx in oracle_case_indices("TY_tables.json", data)
        tab = json_int_matrix(data["Input"][idx])

        check_mt_equal(
          TY_fusion_ring(tab),
          data["Output"][idx],
          "TY_fusion_ring did not match Anyonica TY_tables.json case $idx",
        )
      end
    end

    if want_level("basic") || want_level("intermediate")
      @testset "invalid input" begin
        check_throws(
          () -> TY_fusion_ring(bad_tab), "TY_fusion_ring accepted an invalid group table"
        )
      end
    end
  end

  # ============================================================
  # HI_fusion_ring
  # ============================================================

  @testset "HI_fusion_ring" begin
    z2_tab = [
      1 2
      2 1
    ]

    #= z3_tab = [
      1 2 3
      2 3 1
      3 1 2
    ]
    =#

    D3_tab = [
      1 2 3 4 5 6 
      2 3 1 6 4 5
      3 1 2 5 6 4
      4 5 6 1 2 3
      5 6 4 3 1 2
      6 4 5 2 3 1
    ]

    nonsymmetric_group_tab = D3_tab
    bad_tab = [
      1 1
      2 2
    ]

    maybe_testset("basic", "1. basic construction") do
      r = HI_fusion_ring(z2_tab)

      check_true(r isa FusionRing, "HI_fusion_ring(Z2 table) did not return a FusionRing")
      check_equal(
        size(multiplication_table(r)),
        (rank(r), rank(r), rank(r)),
        "HI_fusion_ring(Z2) did not produce a cubic multiplication table",
      )
      return check_equal(
        length(labels(r)),
        rank(r),
        "HI_fusion_ring(Z2) did not produce the correct number of labels",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      r = HI_fusion_ring(z2_tab)

      check_equal(
        rank(r), 4, "HI_fusion_ring on the Z2 Cayley table did not produce rank 4"
      )
      check_equal(
        size(multiplication_table(r)),
        (4, 4, 4),
        "HI_fusion_ring(Z2) multiplication table had wrong size",
      )
      check_equal(
        labels(r), ["1", "2", "ρ₁", "ρ₂"], "HI_fusion_ring(Z2) labels were incorrect"
      )

      mt1 = multiplication_table(r)

      check_equal(mt1[1, 1, 1], 1, "vacuum × vacuum was not vacuum in HI_fusion_ring(Z2)")
      return check_equal(mt1[1, 3, 3], 1, "vacuum × ρ_1 was not ρ_1 in HI_fusion_ring(Z2)")
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("HI_tables.json")

      for idx in oracle_case_indices("HI_tables.json", data)
        tab = json_int_matrix(data["Input"][idx])

        check_mt_equal(
          HI_fusion_ring(tab),
          data["Output"][idx],
          "HI_fusion_ring did not match Anyonica HI_tables.json case $idx",
        )
      end
    end

    if want_level("basic") || want_level("intermediate")
      @testset "invalid input" begin
        check_throws(
          () -> HI_fusion_ring(bad_tab),
          "HI_fusion_ring accepted a table that was not a valid group table",
        )

        check_throws(
          () -> HI_fusion_ring(nonsymmetric_group_tab),
          "HI_fusion_ring accepted a valid group table that was not symmetric as a matrix",
        )
      end
    end
  end
end
