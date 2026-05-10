# ROGUI Developer Guide

ROGUI is the Resonant Orbit GUI for exploring the nested search space used by this repository's resonant-arc model. It is a MATLAB developer tool, not a separate application: it calls the same patched-conics, Lambert, gravity-assist, and DSM routines used by the scripts in the repository root.

## Launch

From the repository root:

```matlab
addpath(genpath(pwd))
TrajectorySearchSpaceGUI
```

The default case mirrors the Galileo-style `Earth,Venus,Earth,Earth,Jupiter` scenario used in `TEST_MGA_Plot.m`.

## What The GUI Explores

ROGUI separates the resonant optimization into two visible layers:

- Outer map: outgoing hyperbolic excess velocity after the resonant flyby, parameterized by `|vinf_out|`, `theta`, and `phi`.
- Inner map: DSM anomaly `nu_DSM` and branch choices for `revs_before`, `lw`, and `lp`.

Each outer point runs an inner optimization unless `lw` or `lp` are fixed in the GUI. The outer map color is the optimized local resonant-leg cost. Selecting an outer point recomputes the branch-resolved inner map for that fixed `vinf_out`. Selecting an inner point reconstructs the associated trajectory.

## User-Facing Units

The GUI displays angle inputs and plot axes in degrees. Internally, all numerical helpers use radians to remain consistent with the core MATLAB code.

Distances are km in computations and AU in trajectory plots. Velocities are km/s. Public mission TOFs are days; internal Lambert and DSM timing uses seconds.

## Folder Structure

```text
ROGUI/TrajectorySearchSpaceGUI.m
```

Main programmatic `uifigure` entry point. Owns UI layout, callbacks, GUI state, click selection, arrow-key navigation, caching, and metrics refresh.

```text
ROGUI/config/
```

Builds a validated mission/search config. `createSearchSpaceGuiConfig.m` computes planet states, detects resonant legs, prepares adjacent Lambert arcs, and derives the `vinf_i/vinf_o` context for the selected resonant transfer.

```text
ROGUI/search/
```

Builds search-map tables. `runOuterSearchSpaceMap.m` sweeps `|vinf_out|/theta/phi`; `runInnerNuSearchMap.m` evaluates branch-resolved `nu_DSM` samples for a selected outer point.

```text
ROGUI/evaluation/
```

Prepares DSM state models and reconstructs selected trajectories. `evaluateTrajectoryCandidate.m` returns a struct compatible with plotting and metric inspection.

```text
ROGUI/plot/
```

Reusable axes-based plot helpers for the outer map, inner map, and selected trajectory. They contain no optimizer logic.

```text
ROGUI/cache/
ROGUI/validation/
```

Cache-key generation and lightweight smoke checks.

## Typical Development Workflow

1. Run `runGuiSmokeChecks()` after changing wrappers or plotting helpers.
2. Open `TrajectorySearchSpaceGUI`.
3. Keep grids coarse while testing behavior.
4. Use `lw branch` and `lp branch` as `Auto` for baseline behavior; set them to `0` or `1` to isolate branch effects.
5. Click or use arrow keys on the outer map to inspect neighboring solutions.
6. Click branch samples in the inner map to force a specific DSM branch into the trajectory view.

## Important Implementation Notes

- GUI state is stored in `fig.UserData` as an `app` struct.
- The trajectory axes preserve user camera orientation between selection updates.
- The trajectory plot visually scales `z` by `10x`; this is display-only.
- Inner-map samples are not connected by lines because multiple branches can produce different costs at the same `nu_DSM`.
- Generated graph exports are opt-in; the GUI does not write into `1-graphs/` automatically.

## Validation

Run:

```matlab
addpath(genpath(pwd))
results = runGuiSmokeChecks()
```

This builds a small default scenario, evaluates a tiny outer map, evaluates an inner map for the best feasible candidate, and attempts trajectory reconstruction.

## Related Documentation

- `README.md`: repository-level overview, scientific model summary, and workflow.
