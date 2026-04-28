# holobiont-popgen

Forward-time individual-based simulation of holobiont evolution: a host
population in which each individual carries its own microbial community,
and where host fitness depends on a weighted contributed of the host's own
genetically-encoded trait and the mean trait of its microbiome.


## Two sets of simulations

The repo contains the code for two related simulation setups. Each is within
its own subfolder, with its own source code, parameter sweep,
and analysis scripts. 

| Folder                  | What it sweeps                                  | Generations | Tasks |
|-------------------------|-------------------------------------------------|------------:|------:|
| `standard/`             | Vertical inheritance × env-change × *fixed* initial importance values (30 levels from 0 to 0.5). Importance does not mutate. | 3000        | 240   |
| `importance_mutation/`  | Vertical inheritance × env-change × selective vs neutral, with importance free to mutate from an initial value of 0. | 8000        | 720   |

The two setups are largely identicaly but the `importance_mutation/`
variant adds a `metabolic_cost` argument to the host-fitness function
and a `step_size` parameter to the trait mutator (the metabolic cost is not used).

## Repository layout

```
holobiont-popgen/
├── README.md
├── LICENSE
├── .gitignore
├── data/
│   ├── README.md
│   └── env_list_rerun.RDS         # bundled env trajectories (3 scenarios)
├── analysis/
│   ├── analysis_report.Rmd        # ★ knit to reproduce the published figures
│   └── test_simulation.Rmd        # ★ knit to run a small live simulation
├── standard/                      # standard sweep (verbatim code)
│   ├── corefunctions.R
│   ├── *.cpp                      # Rcpp kernels
│   ├── holobiont_setup_array.R
│   ├── holobiont_simulation_array.R
│   ├── holobiont_array.sl
│   ├── process_simresults.R
│   ├── process_simresults.sl
│   ├── diversity_env_extract_array*.sl
│   ├── rewritten_analysiscode_Feb10.R
│   └── suppfigs.R
└── importance_mutation/           # importance-mutation sweep (verbatim code)
    ├── corefunctions.R
    ├── *.cpp                      # Rcpp kernels
    ├── holobiont_setup_array.R
    ├── holobiont_simulation_array.R
    ├── holobiont_array.sl
    ├── analysis_holobiont.sl
    ├── effective_sel_strength.sl
    ├── effective_sel_plots.R
    └── plots_of_importance_evolution.R
```

## What's bundled

Only one binary file is bundled in the repo: `data/env_list_rerun.RDS`,
a 1.2 MB named list of three environment-trajectory matrices (`step_gen1`,
`step_gen5`, `step_gen50`) used by both experiments. The simulation
*outputs* (one ~1 GB RDS per task × ~7,000 tasks total) are *not*
included — they are produced by running the simulations

## Quick start

### Run a small live simulation (no HPC required)

```r
rmarkdown::render("analysis/test_simulation.Rmd")
```

Or open it in RStudio and click **Knit with Parameters…** to tweak
population size, number of generations, etc. With the defaults
(30 hosts × 50 species × 100 generations) it finishes in well under a
minute on a laptop. The Rmd sources the C++/R engine directly from
`importance_mutation/` and produces a small set of diagnostic plots.

### Reproduce the published-figure analysis

```r
rmarkdown::render(
  "analysis/analysis_report.Rmd",
  params = list(
    importance_results_dir   = "/path/to/results_importance_mutation_standard_rerun",
    standard_aggregated_data = "/path/to/combined_data_all_aggregated.RDS"
  )
)
```

The two `params` are paths to the *intermediate* datasets that the
SLURM extraction jobs and the `process_simresults.R` script produce
(see below). If you don't have the standard-experiment dataset,
leave `standard_aggregated_data = ""` and that section is skipped.

### Run the full sweep (HPC)

The two experiment folders each contain a complete pipeline:

1. **Build the parameter registry and job mapping**
   ```bash
   cd standard/                     # or importance_mutation/
   Rscript holobiont_setup_array.R
   ```

2. **Submit the simulation array**
   ```bash
   sbatch holobiont_array.sl        # adjust --account, --array range, paths
   ```
   Each task writes one RDS per (combo × replicate) pair.

3. **Extract per-combo summaries** (only `importance_mutation/`)
   ```bash
   sbatch analysis_holobiont.sl     # variance ratio
   sbatch effective_sel_strength.sl # effective selection
   ```

4. **Aggregate the standard sweep** (only `standard/`)
   ```bash
   sbatch process_simresults.sl     # produces combined_data_all_aggregated.RDS
   ```

5. **Plot**
   - `analysis/analysis_report.Rmd` (recommended), or
   - the original scripts directly:
     `standard/rewritten_analysiscode_Feb10.R`,
     `standard/suppfigs.R`,
     `importance_mutation/effective_sel_plots.R`,
     `importance_mutation/plots_of_importance_evolution.R`.

## Path conventions in the original scripts

The scripts in `standard/` and `importance_mutation/` are **preserved
verbatim** from the project as it was run on NeSI, with absolute
NeSI paths still embedded (e.g. `/nesi/nobackup/uoa04039/...`). On a
fresh clone you'll need to either:

- replace those `setwd(...)` and `Rcpp::sourceCpp(...)` calls with paths
  on your own system, or
- run the scripts from inside the experiment folder with all source
  files alongside (the SLURM scripts handle this with `cd <dir>`).

The two Rmds in `analysis/` are written to be path-portable and don't
require any edits.

## Citation

If you use this code, please cite the accompanying paper (citation TBC)
and link to this repository.

## License

MIT — see [LICENSE](LICENSE).
