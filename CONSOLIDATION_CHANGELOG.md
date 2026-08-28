## Summary

This PR cleans up and reconciles a number of issues across validation, metadata handling, imports, gradings, categorifiability checks, and the test suite.

All changes are marked with `#changed:` to make them easier to review.

## Main changes

* Tightened validation for multiplication tables, inverses, associativity, and subrings.
* Fixed restriction and permutation so labels and other metadata tied to simple objects are preserved correctly.
* Made `force_compute` behavior consistent across dimensions, characters, formal codegrees, decompositions, gradings, and categorifiability checks.
* Fixed several inconsistencies involving tensor product decompositions, upper central series data, universal/all gradings, imports, and automorphism metadata.

## Data and database handling

* Added stricter `JSON3`/`StructTypes` decoding through `RingData`, while keeping the legacy dictionary decoder for bulk and in-memory imports.
* Added typed projective-representation import/export using `data[ogid][name][pair_index]`, including validation for malformed inputs.
* Moved the public ring database lookup functions into `src/database.jl`.
* Simplified `src/FusionRings.jl` so it mainly handles imports, includes, and initialization.
* Replaced writes to the source tree for algebraic-number data with a lazy Scratch cache in a user-writable location.
* Added `is_optimized_db()::Bool`.
* Fixed bulk `Dict` imports.
* Updated compatibility metadata for the Julia 1.12 `Pkg` stdlib series.

## Subrings and injections

* Added the Combinatorics dependency and qualified its usage.
* Fixed `which_injection` so failed permutation matches are skipped safely.
* Implemented `injection_form`.
* Replaced the undefined `sub_ring_tables` placeholder with a typed compatibility error that points users to the supported API.

## Tests and CI

* Added support for `FUSIONRINGS_FORCE_COMPUTE` in the test suite.
* Expanded `test/testdata/README.md` with the fixture layout, Git LFS requirements, test levels/groups, representative-case selection, and commands for running the full oracle tests.
* Updated CI to check out LFS content before testing.
* Representative cases now run by default.
* Added a scheduled/manual full Anyonica oracle run with forced computation enabled.



