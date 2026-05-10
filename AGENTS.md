# Repository Guidelines

## Project Structure & Module Organization

This repository contains MATLAB code for multi-gravity-assist and resonant trajectory analysis.

- `MGA2_PGA2.m` is the main trajectory function combining powered gravity assists, Lambert arcs, and resonant VILM legs.
- `TEST_MGA_Plot.m` is the current executable validation/plotting script.
- `Astro/` contains astrodynamics helpers, ephemeris routines, Lambert solvers, and gravity-assist utilities.
- `0-coreResonantTrajectory/`, `0-solverDSM/`, and `0-solverResonantTrajectory/` contain the resonant transfer and DSM optimization logic.
- `0-costFunctions/`, `0-searchAlgorithms/`, and `0-utilities/` hold optimization wrappers, scan/refine helpers, and reusable trajectory utilities.
- `1-graphs/` stores generated PNG figures. Do not update graph assets unless the numerical or plotting behavior intentionally changed.

## Build, Test, and Development Commands

Run commands from the repository root in MATLAB.

```matlab
addpath(genpath(pwd))
TEST_MGA_Plot
```

Adds all folders to the path and runs the Galileo-style computation and plots.

```matlab
addpath(genpath(pwd))
[jd2k,r,v,vd,va,rpga,dvga,dvdsm,vilm_arcs] = MGA2_PGA2(planets,jd2k0,tofs,N,M);
```

Runs the main API directly with caller-provided planets, dates, transfer times, and resonance values.

## Coding Style & Naming Conventions

Use MATLAB `.m` files with one primary function per file, matching the filename exactly. Keep function names descriptive and existing camelCase/PascalCase conventions intact, for example `computeSingleDSMTransferCost`, `findOptimalDSMParameters`, and `GetBodyICF`.

Use 4-space indentation inside functions and control blocks. Prefer explicit input/output comments with units. Keep numerical constants named near first use, and preserve kilometers, seconds, days, and radians in comments where relevant.

## Testing Guidelines

There is no formal test framework. Use `TEST_MGA_Plot.m` as the smoke test after changing trajectory logic, Lambert calls, DSM optimization, or plotting. For new computational routines, add a small deterministic script or section that exercises a nominal case and one invalid geometry/convergence case. Avoid relying only on plots when scalar outputs such as `dvga`, `dvdsm`, or periapsis radius can be checked.

## Commit & Pull Request Guidelines

Recent commit subjects are short and descriptive, often naming the feature or refactor, such as `FuncionsRefactoritzades` or `Create ResonantVILM2.m`. Keep new commits concise, imperative, and focused on one change.

Pull requests should include the trajectory scenario tested, commands run, relevant output changes, and any regenerated files under `1-graphs/`. Link related issues or thesis task notes when available, and include screenshots only when plot appearance changes.

## Agent-Specific Instructions

Do not reorganize numbered directories without a clear migration reason. Preserve existing MATLAB APIs unless the caller updates are included in the same change.
