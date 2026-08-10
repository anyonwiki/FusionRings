# GitLab weekly testing + FuncDeps architecture checks

This bundle adds a GitLab CI setup for FusionRings.jl with three test levels
and a weekly FuncDeps architecture check.



## What the pipeline does

### Test jobs

The pipeline uses the environment variables already supported by
`test/runtests.jl`:

```julia
FUSIONRINGS_TEST_GROUP
FUSIONRINGS_TEST_LEVEL
```

Jobs included:

```text
test:level1
test:level2
test:level3
test:all
```

`test:level1` runs on pushes, merge requests, and scheduled pipelines.

`test:level2`, `test:level3`, and `test:all` run automatically on scheduled
pipelines. They can also be started manually.

Right now, if the individual test files do not yet call `want_level("level1")`,
`want_level("level2")`, or `want_level("level3")`, then these jobs may run the
same tests. 

## FuncDeps architecture check

The `funcdeps:architecture` job:

1. Installs Graphviz.
2. Installs FuncDeps from GitHub:

```julia
Pkg.add(url="https://github.com/Szagha02/FuncDeps.jl")
```

3. Scans `src/`.
4. Generates:

```text
artifacts/funcdeps/moduledeps.dot
artifacts/funcdeps/moduledeps.svg
artifacts/funcdeps/funcdeps_cross_module.dot
artifacts/funcdeps/funcdeps_cross_module.svg
artifacts/funcdeps/funcdeps_full.dot
artifacts/funcdeps/funcdeps_full.svg
artifacts/funcdeps/interactive_graph.html
artifacts/funcdeps/architecture_report.txt
```

5. Fails the CI job if `ci/architecture_policy.toml` is violated.

## First-pass policy

The first policy enforces the requested rule:

```text
structs should not depend on anything else
```

In precise terms:

```text
functions detected in structs.jl may call other functions in structs.jl,
but may not call discovered project functions from other source files/modules.
```

This may initially fail if `structs.jl` calls helpers such as
`is_constant_array` from `general_functions.jl`. That is useful: it shows the
architecture drift the policy is supposed to catch. To fix that, either move the
helper into `structs.jl`, inline it, or relax the policy by allowing that edge.

## Weekly schedule
(just in case we forget how to build the pipeline - though this can be deleted safely as it has no importance on documentation)
In GitLab:

```text
Build → Pipeline schedules → New schedule
```

Example weekly cron:

```text
0 3 * * 0
```

 Sundays at 03:00 in the schedule timezone.

## Local usage

We can run the architecture tools locally after installing FuncDeps:

```julia
using Pkg
Pkg.activate("ci/FuncDepsEnv"; shared=false)
Pkg.add(url="https://github.com/Szagha02/FuncDeps.jl")
```

Then:

```bash
julia --project=ci/FuncDepsEnv ci/run_funcdeps.jl src artifacts/funcdeps ci/architecture_policy.toml
julia --project=ci/FuncDepsEnv ci/check_architecture.jl src ci/architecture_policy.toml artifacts/funcdeps/architecture_report.txt
```
