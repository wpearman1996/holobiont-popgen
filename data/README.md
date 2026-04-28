# Bundled data

## `env_list_rerun.RDS`

A named list of three matrices encoding the environmental-optimum
trajectory used by both experiments. Each entry corresponds to a
different number of microbial generations per host generation:

| Name         | Shape (rows × cols) | Microbial gens per host gen |
|--------------|--------------------:|----------------------------:|
| `step_gen1`  | 3000 × 1            | 1                           |
| `step_gen5`  | 3000 × 5            | 5                           |
| `step_gen50` | 3000 × 50           | 50                          |

- **Rows** index host generations (1..3000).
- **Columns** index microbial generations within a host generation.
- **Values** are environmental optima on the same scale as
  `traitpool_microbes` in the simulation (roughly `[-2.5, 2.5]`).

### Pattern

The "step" pattern: the optimum sits near `0` (with small noise) for
generations 1–2000, then jumps to ~`0.8` for generations 2001–3000.
This sets up the pre/post-step contrast used throughout the analysis.

The within-row (across-microbial-generation) values interpolate linearly
between successive host-generation optima, so the more microbial
generations there are per host generation, the more gradually the
microbes experience the change within a host's lifetime.

### Origin

Subset of a larger six-element file produced upstream of this repo;
only the first three (`step_gen1`, `step_gen5`, `step_gen50`) are used
by any script in either experiment.
