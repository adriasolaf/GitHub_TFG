# Resonant Orbit Trajectory Toolkit

This repository contains MATLAB code for interplanetary trajectory analysis with multiple gravity assists and optional resonant arcs. The main scientific focus is the resonant `N:M` arc: a same-planet encounter pair where the spacecraft leaves a body, follows heliocentric two-body motion, performs one DSM, and returns after a prescribed planetary phasing interval.

The model is based on patched conics, Keplerian two-body propagation, Lambert transfers, powered gravity-assist estimates, and nested search over resonant-arc variables.

## Repository Layout

```text
MGA2_PGA2.m
TEST_MGA_Plot.m
Astro/
0-coreResonantTrajectory/
0-solverDSM/
0-solverResonantTrajectory/
0-costFunctions/
0-searchAlgorithms/
0-utilities/
ROGUI/
1-graphs/
```

- `MGA2_PGA2.m`: main trajectory API for the patched-conics MGA/VILM workflow.
- `TEST_MGA_Plot.m`: baseline Galileo-style computation and plotting script.
- `Astro/`: astrodynamics helpers, body states, Lambert solver, element/state conversions, and flyby utilities.
- `0-coreResonantTrajectory/`: single DSM transfer construction and cost logic.
- `0-solverDSM/`: inner DSM optimizer over anomaly, revolution split, and Lambert branch.
- `0-solverResonantTrajectory/`: outer resonant optimizer over outgoing `v_infinity`.
- `0-costFunctions/`: resonant total-cost wrappers.
- `0-searchAlgorithms/`: grid scans, scan/refine, and ternary search helpers.
- `0-utilities/`: trajectory reconstruction utilities.
- `ROGUI/`: Resonant Orbit GUI and its script-callable support layer.
- `1-graphs/`: generated plots and figures.

## Quick Start

Run the baseline script:

```matlab
addpath(genpath(pwd))
TEST_MGA_Plot
```

Run the main API directly:

```matlab
planets = {'Earth', 'Venus', 'Earth', 'Earth', 'Jupiter'};
jd2k0 = -3727;
tofs = [115, 301, 730, 1094];
N = 1;
M = 2;

[jd2k,r,v,vd,va,rpga,dvga,dvdsm,vilm_arcs] = MGA2_PGA2(planets,jd2k0,tofs,N,M);
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
- A resonant VILM leg is detected when consecutive bodies in `planets` are equal, for example `Earth -> Earth`.
- Each resonant leg uses one DSM and returns to the same planet after `M` planet revolutions while the spacecraft completes `N` heliocentric revolutions.

For a candidate resonant outgoing velocity:

```text
r_sc0 = r_planet(jd0)
v_sc0 = v_planet(jd0) + vinf_out
tof_total = M * T_p
```

The code converts this state to Keplerian elements, rejects hyperbolic initial heliocentric orbits, propagates to a DSM anomaly, then solves a Lambert arc from the DSM point to the planet position at `jd0 + M*T_p`.

The resonant objective is:

```text
J = DeltaV_GA1 + DeltaV_DSM + DeltaV_GA2
```

where edge resonances omit the missing launch-side or arrival-side flyby term. The optimizer is nested: the outer search varies `vinf_out` magnitude and direction, while the inner search varies `nu_DSM`, revolution split, and Lambert branch flags `lw/lp`.

This is a deterministic patched-conics design model, not a high-fidelity propagator. It excludes third-body perturbations, finite burns, atmosphere, oblateness, solar radiation pressure, and ephemeris uncertainty. The search strategy is heuristic and does not prove a global optimum over all mission sequences, epochs, or TOFs.

## ROGUI Workflow

ROGUI is a developer tool for exploring the resonant search topology:

1. Configure the mission, resonance, and search grids.
2. Generate the outer PCP-like map over `|vinf_out|`, `theta`, and `phi`.
3. Select an outer point by click or arrow keys.
4. Inspect the branch-resolved inner `nu_DSM` map.
5. Select a DSM sample and inspect the resulting trajectory and metrics.

See `ROGUI/README.md` for the developer tutorial and implementation guide.

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

- `AGENTS.md`: contributor guidelines.
- `ROGUI/README.md`: ROGUI developer guide.
