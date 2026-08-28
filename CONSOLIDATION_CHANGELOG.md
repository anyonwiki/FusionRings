# FusionRings.jl Consolidation Changelog

Date: 2026-08-28

Branch: `consolidated_changes`

Base: `develop` (`7f37757200949a5c7a0c32f876ea8ea194531369`)

## Purpose

This branch consolidates the compatible intent of `bugfix_suite` and
`categorifiability_changes` against the current `FusionRing` schema. It also
replays all six dependency-update branches. The existing `main` and `develop`
refs remain unchanged.

Every source function added or modified by this consolidation now has an
immediately preceding `#changed:` comment summarizing the local correction.

## Functional changes

| Area | Summary |
| --- | --- |
| Constructor validation | Reject rank-zero tables; validate a unique matching left/right dual for every simple; validate cached character, FP-dimension, subring, grading, and upper-central-series metadata; fix vector `texnames` handling. |
| Subrings | Reject duplicate subring indices, check closure against the parent table, and preserve selected labels in computed restricted rings. |
| Permutations | Reindex labels, character columns, FP dimensions, subring/UCS indices, grading values, and automorphism data consistently; validate permutation shape and unit preservation. |
| FP dimensions and caching | Select exact/numeric Perron eigenvalues robustly; return real numeric FP values; propagate `force_compute` through scalar/component FP dimensions, sorting, barcode, characters, formal codegrees, decompositions, gradings, and nilpotence. |
| Decompositions | Normalize the supported tensor-product kind, preserve positional compatibility, propagate forced subring discovery, fix the `d1`/`d2` comparison typo, and compare factor multisets by fusion-ring equivalence rather than UUID identity alone. |
| Gradings and central series | Validate universal-grading coverage and unique output grades, restore stored `all_gradings` reuse, and ensure forced recomputation bypasses nested adjoint/UCS caches. |
| Import/export | Validate current JSON shapes and types, handle absent automorphism metadata, normalize the legacy `TensorProduct` key, preserve character-table orientation on round trip, and synchronize the QQBar disk/in-memory cache. |
| Error behavior | Replace undefined `message(...)` paths with explicit `ArgumentError`s for invalid sort choices and non-group Cayley-table requests. |
| Categorifiability | Enable and export CSP, PDC, d-number, zero-spectrum, and one-spectrum criteria; correct CSP indexing and denominator; add exact algebraic-integrality support; repair the d-number divisibility direction; replace undefined ZSC/OSC helpers with project APIs. |
| Tests | Add regression coverage for validation, label/metadata permutations, stale-cache bypass, decompositions, import/export orientation and validation, d-numbers, and categorifiability criteria, including the existing Anyonica reference-test structure. |
| Dependencies | Replay checkout 7, codecov 7, cache 3, setup-julia 3, JSON 1.7.1, and Documenter 1.17.0. |

## Where Ostrik d-numbers enter

Ostrik proves that every formal codegree of a (multi-)fusion category is a
d-number, although the formal codegrees of a general based ring need not be.
Therefore, a fusion ring with even one formal codegree that is not a d-number
cannot be the Grothendieck ring of a fusion category. This is a necessary
condition, not a sufficiency test: passing it does not prove categorifiability.

The implementation path is:

1. `dn_criterion(fr; force_compute=...)` obtains `formal_codegrees(fr)`.
2. Each formal codegree is passed to `is_d_number`.
3. `is_d_number` first requires algebraic integrality.
4. For the monic minimal polynomial
   `x^n + a_1*x^(n-1) + ... + a_n`, it checks Ostrik's equivalent coefficient
   condition that `a_i^n / a_n^i` is an integer.
5. If any formal codegree fails, `dn_criterion` returns `true`, consistently
   with the package convention that `true` means “the criterion found an
   obstruction.”

The current executable criterion is restricted to commutative rings. Ostrik's
theorem itself is not commutativity-restricted; the implementation guard exists
because the current noncommutative matrix path yields the central-element
eigenvalue `f_E * dim(E)`, and does not yet recover `dim(E)` well enough to
extract every `f_E` reliably.

Reference: Victor Ostrik, [*On formal codegrees of fusion categories*](https://arxiv.org/abs/0810.3242), especially Theorem 1.2 and Lemma 2.7.

## Compatibility and semantics

- No exported source function was intentionally deleted.
- The internal Boolean helpers `crit1`, `crit2`, and `crit3` were replaced by
  raw-sum helpers plus explicit Boolean wrappers so ZSC and OSC can test the
  distinct zero/one conditions correctly.
- The positional `fpdim(ring, true)` form remains supported.
- Positional decomposition-kind calls remain supported.
- Categorifiability criteria consistently use obstruction semantics:
  `true` means the criterion rules out the relevant categorification; `false`
  means only that this criterion found no obstruction.
- Malformed stored JSON now fails early instead of being silently normalized
  into inconsistent metadata.

## Commit groups

| Commit | Contents |
| --- | --- |
| `00b9fe9` | Validation, duality, and subring invariants |
| `c8d0348` | Permutations, labels, and cached metadata |
| `0bf51dd` | Exact/numeric FP dimensions and forced recomputation |
| `a2f13c9` | Decompositions, gradings, imports, TeX names, and QQBar persistence |
| `ba3ed05` | Categorifiability criteria and tests |
| `a01e98d`–`e11ab8f` | Six dependency updates |
| `77d36df` | CI line-ending normalization |
| `ed267a2` | Nested cache propagation and final metadata checks |
| Current tip | Per-function `#changed:` annotations and this professor-facing changelog |

## Verification status

Static verification includes clean Git diffs, TOML/YAML parsing, conflict-marker
and unresolved-helper scans, UTF-8/NUL checks, bundle integrity, ancestry checks,
and a round-trip fetch into a fresh repository.

Julia runtime tests remain the final approval gate because the consolidation
workspace does not contain a Julia executable or the Git LFS payloads for the
Anyonica fixtures. In a full clone, run `git lfs pull` followed by
`julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.test()"` before merging
into `develop`.
