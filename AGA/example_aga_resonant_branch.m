%% Example AGA for a resonant branch
% This script mirrors example_aga_rastrigin.m but uses the resonant cost.
%
% AGA individual format:
%   [vinfo_km_s, theta_deg, phi_deg, nu_deg, lp]
%
% The chromosome controls the outgoing excess velocity, the DSM anomaly
% sample, and the Lambert long-period flag. The cost function still selects
% the best feasible revs_before internally because the current resonant
% solver requires that discrete split.

%% Scenario setup

% Use the same Galileo-style sequence used by the repository smoke tests.
% The selected resonant branch is Earth -> Earth at index 3.
planets = {'Earth', 'Venus', 'Earth', 'Earth', 'Jupiter'};
jd2k0 = -3727;
Nres = 1;
Mres = 2;
resonantIndex = 3;

% Store only the fields required by aga_fun_cost. Distances are km,
% velocities are km/s, epochs are JD2000 days, and periods are seconds.
scn.mu_sun = GetBodyProps('Sun');
scn.body = planets{resonantIndex};
scn.jd0 = jd2k0 + 115 + 301;
scn.N = Nres;
scn.M = Mres;
scn.apsis_flag = 1;
scn.res_flag = 2;
scn.penaltyCost = 1e12;

% Planet constants and resonant body orbital period.
[scn.mu_planet, bodyRadius] = GetBodyProps(scn.body);
scn.vmr_safety = 1.05 * bodyRadius;
[smaPlanet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(scn.body, scn.jd0);
scn.T_p = 2*pi*sqrt(smaPlanet^3 / scn.mu_sun);

% The resonant leg duration is M planetary revolutions.
tofs = [115, 301, scn.T_p*Mres/86400, 1094];
days2secs = 86400;
jd2k = zeros(numel(planets), 1);
r = zeros(numel(planets), 3);
v = zeros(numel(planets), 3);
vd = NaN(numel(planets) - 1, 3);
va = NaN(numel(planets) - 1, 3);

% Compute encounter epochs and planet states.
jd2k(1) = jd2k0;
for i = 2:numel(planets)
    jd2k(i) = jd2k(i - 1) + tofs(i - 1);
end
for i = 1:numel(planets)
    [r(i, :), v(i, :)] = GetBodyICF(planets{i}, jd2k(i), scn.mu_sun, 1);
end

% Solve only the adjacent non-resonant Lambert arcs needed by the local
% resonant cost: Venus -> Earth before resonance and Earth -> Jupiter after.
for i = [resonantIndex - 1, resonantIndex + 1]
    dnu = DeltaNu3(r(i, :), r(i + 1, :), 1);
    lw = double(dnu > pi);
    [vd(i, :), va(i, :)] = Lambert(r(i, :), r(i + 1, :), tofs(i)*days2secs, scn.mu_sun, lw, 0, 0);
end

% Store the incoming and required outgoing excess velocities for the
% resonant branch. These define the two powered flyby matching terms.
scn.vinfi_req = va(resonantIndex - 1, :) - v(resonantIndex, :);
scn.vinfo_req = vd(resonantIndex + 1, :) - v(resonantIndex + 1, :);

% Search ranges (min, max). 
% Rows match the GA chromosome format.
ranges = [
    0.5, 8.0      % vinfo [km/s]
   -180, 180      % theta [deg]
    0, 0          % phi [deg]
   -180, 180      % nu [deg]
    0, 1          % lp [-]
];

% Mutation standard deviations (sigma) for [vinfo, theta, phi, nu].
mutationStep = [0.7, 20, 0, 20];

%% AGA setup

% AGA options follow the same convention as example_aga_rastrigin.m.
opts.ninfo = 2;
opts.label = 20;
opts.dopar = 0;
opts.nhist = 1;

% Target cost (m/s)
goal = 1e-3; 

% Population configuration
ng = 15; % Number of generations
np = 40; % Population size
Npop = [
    2, ...              % Elites
    floor(np*0.20), ... % Mutants
    floor(np*0.10), ... % Newcomers
    floor(np*0.30)      % Parent Pool Size
];

% AGA callback function definitions. 
% Keeping these inline makes it clear which operator each AGA slot uses.
unifun = @(pop, fit) deal(pop, fit);
fitfun = @(ind) aga_fun_cost(ind, scn, ranges);
mutfun = @(ind, fit) aga_fun_mutate(ind, fit, ranges, mutationStep);
repfun = @(a, b, fa, fb) aga_fun_crossover(a, b, fa, fb, ranges);
ranfun = @() aga_fun_random(ranges);
prifun = @(ind) aga_fun_print(ind);

%% Run

% The heuristic is intentionally non-repeatable unless a fixed seed is set.
rng('shuffle');

% Call AGA optimizer
[bestInd, bestFit, nite, lastPop, lastFit, history] = aga( ...
    opts, np, goal, ng, Npop, unifun, fitfun, ...
    mutfun, repfun, ranfun, prifun);

% Print best individual found
fprintf('\nBest resonant AGA individual after %d generations:\n', nite);
aga_fun_print(bestInd);
fprintf('\nFitness = %.6f km/s\n', bestFit);

% Re-evaluate the winner to inspect the debug breakdown.
[~, dbg] = aga_fun_cost(bestInd, scn, ranges);
disp(dbg);
