# Changelog (just to see what changed)


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



## Some bugfixes:

- No exported source function was intentionally deleted.
-  `crit1`, `crit2`, and `crit3` were replaced by
  raw-sum helpers plus explicit Boolean wrappers so ZSC and OSC can test 
  distinct zero/one conditions correctly.
- Positional decomposition-kind calls remain supported.
- Malformed stored JSON now fails early instead of being silently normalized
  into inconsistent metadata.

## Commit groups

I grouped the commits into several categories so hopefully that allow to see changes nicely

| Commit | Contents |
| --- | --- |
| Validation, duality, and subring invariants |
| Permutations, labels, and cached metadata |
|  Exact/numeric FP dimensions and forced recomputation |
| Decompositions, gradings, imports, TeX names, and QQBar persistence |
| Categorifiability criteria and tests |
|  Six dependency updates |
| CI line-ending normalization |
|  Nested cache propagation and final metadata checks |



