function results = runGuiSmokeChecks()
%runGuiSmokeChecks ROGUI helper smoke checks
%   Runs a small script-callable validation of config creation, outer map
%   evaluation, inner map evaluation, and trajectory reconstruction.
%
% Inputs:
%   [-]
%
% Outputs:
%   results: validation result structure
%
% Example:
%   [ results ] = runGuiSmokeChecks ();
%
% References:
%   [-]
%
%May 2026

    addpath(genpath(pwd));

    results = struct();
    results.ok = false;
    results.messages = strings(0, 1);

    try
        % Keep the grid tiny so this check can be run frequently.
        cfgInput = struct();
        cfgInput.outerVmag = linspace(2.0, 4.0, 2);
        cfgInput.outerTheta = linspace(-0.5, 0.5, 2);
        cfgInput.outerPhi = 0;
        cfgInput.innerNu = linspace(0.05, pi - 0.05, 8);
        cfgInput.innerScanPoints = 20;
        cfgInput.innerRefineIterations = 10;

        config = createSearchSpaceGuiConfig(cfgInput);
        results.messages(end + 1) = "Config created.";

        [outerMap, bestOuter] = runOuterSearchSpaceMap(config);
        assert(~isempty(outerMap), 'Outer map is empty.');
        results.messages(end + 1) = "Outer map evaluated.";

        if isempty(bestOuter)
            results.messages(end + 1) = "No feasible outer candidate in tiny smoke grid.";
            results.ok = true;
            return;
        end

        [innerMap, bestInner] = runInnerNuSearchMap(config, bestOuter);
        assert(~isempty(innerMap), 'Inner map is empty.');
        results.messages(end + 1) = "Inner map evaluated.";

        if isempty(bestInner)
            traj = evaluateTrajectoryCandidate(config, bestOuter);
        else
            traj = evaluateTrajectoryCandidate(config, bestOuter, bestInner);
        end
        assert(isfield(traj, 'status'), 'Trajectory result has no status field.');
        results.messages(end + 1) = "Trajectory candidate evaluated.";

        results.ok = true;
    catch err
        results.ok = false;
        results.messages(end + 1) = string(err.message);
    end
end
