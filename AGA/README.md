# AGA Developer Guide

## Purpose

`AGA/` contains a simple in-house Genetic Algorithm optimizer and the repository-specific operators used to apply it to resonant orbit branch optimization.

The generic AGA engine minimizes a scalar fitness function. It does not know about trajectory physics, bounds, integer variables, or invalid geometry. Those rules are supplied through callback functions.

## Generic AGA Logic

Main entry point:

```matlab
[bestind, bestfit, nite, lastpop, lastfit, history] = aga( ...
    opts, pop, goal, ng, N, unifun, fitfun, ...
    mutfun, repfun, ranfun, prifun);
```

AGA minimizes `fitfun(ind)`.

Each generation:

1. Calls `unifun` to remove repeated individuals if desired.
2. Refills missing population slots with `ranfun`.
3. Evaluates every individual with `fitfun`.
4. Sorts the population by increasing fitness.
5. Saves history according to `opts.nhist`.
6. Stops if `bestfit <= goal` or the generation limit is reached.
7. Builds the next generation as elites, mutants, descendants, and newcomers.

Population control is:

```matlab
N = [ne, nm, nn, na];
```

- `ne`: elite individuals copied unchanged.
- `nm`: mutants created with `mutfun`.
- `nn`: random newcomers created with `ranfun`.
- `na`: number of best individuals eligible as parents.
- descendants fill the remaining population slots using `repfun`.

## Generic Callback Contract

AGA expects these callbacks:

```matlab
unifun(pop, fit) -> [pop, fit]
fitfun(ind) -> scalar
mutfun(ind, fit) -> ind
repfun(parentA, parentB, fitA, fitB) -> ind
ranfun() -> ind
prifun(ind)
```

Individuals are stored in a cell array. The individual itself can be any MATLAB object, but this repository uses numeric vectors.

## Resonant Orbit AGA Implementation

The resonant branch chromosome is:

```matlab
ind = [vinfo_km_s, theta_deg, phi_deg, nu_deg, lp];
```

Gene meanings:

- `vinfo_km_s`: outgoing hyperbolic excess speed [km/s].
- `theta_deg`: outgoing excess velocity azimuth [deg].
- `phi_deg`: outgoing excess velocity elevation [deg].
- `nu_deg`: DSM anomaly parameter [deg].
- `lp`: Lambert long-period branch flag, `0` or `1`.

Default broad ranges are:

```matlab
ranges = [
    0, 20
   -180, 180
   -180, 180
   -180, 180
    0, 1
];
```

The current resonant solver still requires `revs_before`. The AGA chromosome does not include it. `aga_fun_cost` keeps `nu` and `lp` fixed from the chromosome, evaluates all admissible `revs_before` values, and uses the best feasible one.

The cost is:

```text
J = dV_GA1 + dV_DSM + dV_GA2
```

Invalid candidates return a large penalty cost instead of throwing into AGA.

## Resonant Scenario Structure

`aga_fun_cost` receives a scenario structure named `scn`.

Required or auto-fillable fields:

- `body`: resonant body name.
- `jd0`: resonant branch departure epoch [JD2000 days].
- `N`: spacecraft revolutions in the resonant orbit.
- `M`: planet revolutions in the resonant orbit.
- `apsis_flag`: DSM reference apsis flag.
- `res_flag`: resonance position flag, `1` launch edge, `2` middle, `3` arrival edge.
- `vinfi_req`: incoming excess velocity before the resonant branch [km/s].
- `vinfo_req`: required outgoing excess velocity after the resonant branch [km/s].
- `mu_sun`: Sun gravitational parameter [km^3/s^2].
- `mu_planet`: resonant body gravitational parameter [km^3/s^2].
- `vmr_safety`: minimum flyby radius [km].
- `T_p`: resonant body orbital period [s].
- `penaltyCost`: optional invalid-candidate penalty.

`example_aga_resonant_branch.m` shows how to build this structure directly for the Galileo-style `Earth -> Earth` branch.

## Typical Resonant AGA Usage

```matlab
addpath(genpath(pwd))

ranges = [0.5 8; -180 180; 0 0; -180 180; 0 1];
mutationStep = [0.7, 20, 0, 20];

unifun = @(pop, fit) deal(pop, fit);
fitfun = @(ind) aga_fun_cost(ind, scn, ranges);
mutfun = @(ind, fit) aga_fun_mutate(ind, fit, ranges, mutationStep);
repfun = @(a, b, fa, fb) aga_fun_crossover(a, b, fa, fb, ranges);
ranfun = @() aga_fun_random(ranges);
prifun = @(ind) aga_fun_print(ind);

[bestInd, bestFit] = aga(opts, np, goal, ng, Npop, ...
    unifun, fitfun, mutfun, repfun, ranfun, prifun);
```

## Folder Contents

Generic optimizer:

- `aga.m`: generic Genetic Algorithm minimizer.
- `example_aga_rastrigin.m`: original numerical example.

Resonant branch example:

- `example_aga_resonant_branch.m`: explicit Galileo-style resonant branch example.

Resonant AGA operators:

- `aga_fun_cost.m`: resonant branch fitness function.
- `aga_fun_random.m`: random individual generator.
- `aga_fun_mutate.m`: mutation operator.
- `aga_fun_crossover.m`: crossover operator.
- `aga_fun_decode.m`: chromosome clamp/decode helper.
- `aga_fun_unique.m`: optional near-duplicate population filter.
- `aga_fun_print.m`: compact chromosome printer.

Shared auxiliary helpers:

- `aga_aux_default_ranges.m`
- `aga_aux_clamp.m`
- `aga_aux_empty_data.m`
- `aga_aux_empty_dbg.m`
- `aga_aux_duplicate_key.m`
- `aga_aux_fitness_at.m`
- `aga_aux_rand_in_range.m`

## Validation

Basic checks from the repository root:

```matlab
addpath(genpath(pwd))
run('AGA/example_aga_resonant_branch.m')
results = runGuiSmokeChecks()
```

`example_aga_resonant_branch.m` is stochastic. The exact best fitness can change between runs unless `rng(seed)` is fixed.
