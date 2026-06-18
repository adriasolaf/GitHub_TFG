# Resonant Orbit Trajectory Toolkit

This repository contains MATLAB code for interplanetary trajectory analysis with multiple gravity assists and optional resonant arcs. The main scientific focus is the resonant `N:M` arc: a same-planet encounter pair where the spacecraft leaves a body, follows heliocentric two-body motion, performs one DSM, and returns after a prescribed planetary phasing interval.

The model is based on patched conics, Keplerian two-body propagation, Lambert transfers, powered gravity-assist estimates, and nested search over resonant-arc variables.

## Repository Layout

```text
MGA2_PGA2.m
TEST_MGA_Plot.m
0-coreResonantTrajectory/
0-solverDSM/
0-solverResonantTrajectory/
0-costFunctions/
0-searchAlgorithms/
0-utilities/
Astro/
AGA/
ROGUI/
```

- `MGA2_PGA2.m`: main trajectory API for the patched-conics MGA/VILM workflow.
- `TEST_MGA_Plot.m`: baseline JUICE-style computation and plotting script.
- `0-coreResonantTrajectory/`: single DSM transfer construction and cost logic.
- `0-solverDSM/`: inner DSM optimizer over anomaly, revolution split, and Lambert branch.
- `0-solverResonantTrajectory/`: outer resonant optimizer over outgoing `v_infinity`.
- `0-costFunctions/`: resonant total-cost wrappers.
- `0-searchAlgorithms/`: grid scans, scan/refine.
- `0-utilities/`: trajectory reconstruction utilities.
- `Astro/`: astrodynamics helpers, body states, Lambert solver, element/state conversions, and flyby utilities.
- `AGA/`: Genetic Algorithm optimizer for resonant-branch search.
- `ROGUI/`: Resonant Orbit GUI and its script-callable support layer.

## Quick Start

Run the baseline script:

```matlab
addpath(genpath(pwd))
TEST_MGA_Plot
```

Launch the GUI:

```matlab
addpath(genpath(pwd))
TrajectorySearchSpaceGUI
```

## Scientific Model

The code uses patched conics and impulsive maneuvers:

- Between encounters, the spacecraft follows heliocentric two-body Keplerian motion under `mu_sun`.
- At encounters, the planet-relative excess velocity is `vinf = v_spacecraft - v_planet`.
- Non-resonant legs are Lambert arcs between planet positions.
- Powered flybys estimate the periapsis impulse needed to turn `vinf_in` into `vinf_out`, with a `1.05 * body_radius` minimum periapsis safety limit.
- A resonant VILM leg is detected when consecutive bodies in `planets` are equal, for example `Earth -> Earth`, and the transfer time matches `M*T_p` within the configured tolerance.
- Each resonant leg uses one DSM and returns to the same planet after `M` planet revolutions while the spacecraft completes `N` heliocentric revolutions.

For a candidate resonant outgoing velocity:

```text
r_sc0 = r_planet(jd0)
v_sc0 = v_planet(jd0) + vinf_out
tof_total = M * T_p
```

The code converts this state to Keplerian elements, rejects hyperbolic initial heliocentric orbits, propagates to a DSM anomaly, then solves a Lambert arc from the DSM point to the planet position at `jd0 + M*T_p`.

The DSM point is parameterized by `nu_DSM`, a revolution split `revs_before`, and the apsis convention:

```text
apsis_flag = 1: nu_m = pi - nu_DSM     % apoapsis reference
apsis_flag = 0: nu_m = -nu_DSM         % periapsis reference
mr_lambert = N - revs_before - 1
```

The canonical inner search scans `nu_DSM` over `[0, 2*pi]`, all `revs_before = N-1 ... 0`, and both Lambert long-period branches when admissible. True-anomaly to mean-anomaly conversion uses a quadrant-safe helper, and flyby turn-angle cosine values are clamped before `acos` to avoid numerical false invalids.

The resonant objective is:

```text
J = DeltaV_GA1 + DeltaV_DSM + DeltaV_GA2
```

where edge resonances omit the missing launch-side or arrival-side flyby term. The optimizer is nested: the outer search varies `vinf_out` magnitude and direction, while the inner search varies `nu_DSM`, revolution split, and Lambert long-period flag `lp`. Lambert transfer direction `lw` is auto-selected internally for the resonant DSM-to-arrival arc.

This is a deterministic patched-conics design model, not a high-fidelity propagator. It excludes third-body perturbations, finite burns, atmosphere, oblateness, solar radiation pressure, and ephemeris uncertainty. The search strategy is heuristic and does not prove a global optimum over all mission sequences, epochs, or TOFs.

## Resonant Model Implementations

The repository has three front ends around the same resonant physical model:

- Main grid/refine optimizer: `MGA2_PGA2.m` with `0-*` solver and cost folders. This is the canonical implementation.
- AGA: a Genetic Algorithm search in `AGA/` that evaluates fixed chromosomes through the same DSM kernel and flyby cost model.
- ROGUI: a GUI and script-callable search-space inspection layer in `ROGUI/`.

Current equivalence status:

| Item | Main grid/refine | AGA | ROGUI |
|---|---:|---:|---:|
| Same-body resonant leg | yes | scenario-provided | yes |
| Checks TOF against `M*T_p` before classifying resonance | yes | scenario-provided | yes |
| Resonant arrival epoch uses `M*T_p` | yes | yes | yes |
| Initial state `v_p0 + vinf_out` | yes | yes | yes |
| Rejects hyperbolic initial heliocentric orbit | yes | yes | yes |
| Uses shared `computeSingleDSMTransferCost` | yes | yes | yes |
| Searches all `revs_before` splits | yes | yes | yes |
| Searches `nu_DSM` over `[0, 2*pi]` | yes | chromosome-dependent, ranges cover full circle | yes |
| Uses cosine clamp for flyby `acos` | yes | yes | yes |
| Uses quadrant-safe true-to-mean anomaly conversion | yes | yes | yes |

The prior ROGUI divergence has been corrected. ROGUI now uses the same same-body plus `M*T_p` resonance classification as the main solver, full-domain DSM anomaly coverage, and the same patched-conics DSM kernel. It still duplicates some scientific glue logic for GUI-specific fixed-branch exploration and trajectory reconstruction; future model changes should either update those paths together or move branch enumeration and flyby-cost helpers into shared utilities.

## ROGUI Workflow

ROGUI is a developer tool for exploring the resonant search topology:

1. Configure the mission, resonance, and search grids.
2. Generate the outer PCP-like map over `|vinf_out|`, `theta`, and `phi`.
3. Select an outer point by click or arrow keys.
4. Inspect the branch-resolved inner `nu_DSM` map.
5. Select a DSM sample and inspect the resulting trajectory and metrics.

See `ROGUI/README.md` for the developer tutorial and implementation guide.

Default GUI grids are:

- `|vinf|`: `1:0.5:10` km/s.
- `theta`: `-180:10:180` deg.
- `phi`: `0` deg.
- `nu_DSM`: `0:3.025:360` deg.

The inner branch plot uses the selected color metric as both its y-value and color scale. The trajectory plot displays positions in AU with equal x/y scale and visual-only `10x` z exaggeration.

## Validation

There is no formal test framework. Use these smoke checks:

```matlab
addpath(genpath(pwd))
TEST_MGA_Plot
```

```matlab
addpath(genpath(pwd))
results = runGuiSmokeChecks()
```

For scientific changes, compare scalar outputs such as total delta-v, `dvga`, `dvdsm`, periapsis radii, selected `nu_DSM`, and selected Lambert branch flags. Do not rely only on plot appearance.

## Development Notes

- Keep GUI-owned code inside `ROGUI/`.
- Preserve existing scientific APIs unless caller updates are included.
- Use MATLAB `.m` files with one primary function per file.
- Distances are km, velocities are km/s, public TOFs are days, and internal Lambert/DSM times are seconds.
- GUI-facing angles are degrees; numerical helpers use radians.

## Documentation

- `ROGUI/README.md`: ROGUI developer guide.
