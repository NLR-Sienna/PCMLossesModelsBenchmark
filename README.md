# PCMLossesModelsBenchmark

Benchmark of transmission-loss and circular-flow modeling approaches built on the
[Sienna](https://www.nrel.gov/analysis/sienna.html) / PowerSimulations.jl stack.

## Examples & documentation

New here? Start with these READMEs:

- **[`Systems/`](Systems/README.md)** — the system builders (5-bus, RTS, CATS) that
  every script consumes.
- **[`SiennaScripts/`](SiennaScripts/README.md)** — the reusable model/simulation
  library (build, solve, post-process).
  - **[`SiennaScripts/CircularFlows/`](SiennaScripts/CircularFlows/README.md)** —
    detect circular (loop) flows, with a catalog of runnable examples. Best
    first run: `julia SiennaScripts/CircularFlows/scripts/rts_example_local.jl`.
  - **[`SiennaScripts/FlowCancelling/`](SiennaScripts/FlowCancelling/README.md)** —
    transmission-investment models with Big-M flow-cancelling terms.
