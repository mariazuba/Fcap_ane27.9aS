# Management Strategy Evaluation for anchovy ane.27.9aS

This repository contains the conditioning, simulation, validation and performance evaluation scripts for European anchovy in the Gulf of Cadiz.

The MSE evaluates management procedures based on an Fcap and biomass escapement rule, including assessment error and forecast uncertainty.

## Simulation design

- Historical period: 1989–2024.
- Projection period: 2025–2054.
- 100 bootstrap-conditioned operating models.
- 10 replicate trajectories per operating model.
- 1000 trajectories per management procedure and configuration.
- 10 simulation blocks of 100 trajectories.
- Four seasons per year.
- Spawning biomass in Q2 and recruitment in Q3.

The robustness design includes REF, R-OM, S-C1, S-C2, I0 and I1. Existing reference simulations are reused where appropriate.

## Repository structure

| Location | Contents |
|----------|----------|
| `boot/data/` | Base SS3 model and source inputs |
| `data/` | Bootstrap inputs, conditioned objects and simulation data |
| `functions/` | HCR and observation implementations |
| `scripts/` | Numbered workflow scripts, diagnostics and validations |
| `cesga/` | Slurm submission scripts |
| `outputs/` | Evaluation tables, figures and diagnostics |
| `report/WD/` | Working Document sources and PDF reports |
| `report/reproducibility/` | Recorded software environments |

## Main workflow files

- `data.R`: preparation of OM and MP inputs.
- `model.R`: scenario definitions and simulation runs.
- `output.R`: performance evaluation and MP selection.
- `report.R`: figures and Working Document rendering.

Calls are commented out by default. Uncomment only the required stages and run them from the project root.

These files do not automatically resolve dependencies. Several original scripts clear the R workspace, so use a fresh R session when executing a stage.

## Workflow order

1. Prepare bootstrap-conditioned inputs using scripts 01–10. Script 02 calls the SS3-to-FLR converter in script 03.
2. Run and validate the reference OM using scripts 12 and 13.
3. Prepare assessment error and forecast uncertainty using scripts 14 and `15_SS3_FORECAST_UNCERTAINTY.R`.
4. Run the relevant observation and HCR validations.
5. Define, simulate and evaluate candidate MPs using scripts 19–24, including the refined-grid stages where required.
6. Select MPs for robustness using script 25.
7. Define robustness configurations using script 26 and run them using script 27.
8. Generate figures and render Working Documents after their required inputs are available.

Diagnostic and validation scripts remain in `scripts/`. Their outputs are also required by several Working Documents.

## Operational functions

`FcapBpaHCR_ane9aS_PROBABILITY_OPTION.R` is the operational HCR used by the candidate MP and robustness runners.

`FcapBpaHCR_ane9aS_CORRECTED.R` is retained for the reference comparison in validation script 18.

`perfectObs4seas.R` supports seasonal observation processing.

## Software requirements

R, the packages required by each script, and the SS3 executable are required. PDF reports additionally require Pandoc and LaTeX.

Recorded local and CESGA package versions are available in `report/reproducibility/`.

The CESGA robustness launcher uses R 4.4.2, one CPU, 16 GB of memory and a four-hour time limit per block.

## CESGA execution

Submit jobs from the project root. Create the log directory first:

```bash
mkdir -p outputs/mse/logs
sbatch cesga/RUN_MP_ROBUSTNESS.sh S-C1 10 full
```

This submits ten blocks for reference MP number 10 under S-C1. Check both Slurm completion status and the saved simulation outputs.

## Data availability and reproducibility

Conditioned objects in `data/Rdata/` provide inputs for running
the MSE.

Full simulation results are excluded from GitHub:

- `data/mse/MP_runs/`
- `data/mse/robustness_runs/`

These results are required to recalculate performance indicators
without rerunning the simulations. They can be generated using
scripts 20 and 27 with the corresponding inputs and configurations.

Input objects and simulation results are planned to be shared
through Zenodo. Download links and extraction instructions will
be added once the dataset is published.

Rebuilding the conditioned objects from the original data requires
the SS3 inputs and intermediate runs described in each script.
Some intermediate runs are not included in this repository.

Recorded R and package versions are available in
`report/reproducibility/`. Scripts that generate random errors
include fixed seeds.

The complete workflow has not yet been verified from a clean
copy of the repository.