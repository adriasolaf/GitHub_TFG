function [outerSel, innerSel] = computeRefinedOptimum(config, bestOuter)
%   Takes the best (vmag, theta, phi) found by runOuterSearchSpaceMap and
%   applies fminsearch directly so the refinement keeps the basin found by
%   the GUI's fine outer grid. 
%
% Inputs:
%   config: normalized ROGUI mission/search configuration structure
%   bestOuter: best row/struct from runOuterSearchSpaceMap (must contain
%              vmag, theta, phi)
%
% Outputs:
%   outerSel: struct with theta, phi, vmag, totalDv, isFeasible
%   innerSel: struct with nu, totalDv, dV_DSM, revs_before, lp, isFeasible
%
% References:
%   [-]
%
% Adria Sola Foixench
% May 2026

    outerSel = [];
    innerSel = [];
    if isempty(bestOuter)
        return;
    end

    % Pull starting point from the GUI's discrete optimum
    if istable(bestOuter)
        bestOuter = bestOuter(1, :);
        vmag0 = bestOuter.vmag;
        theta0 = bestOuter.theta;
        phi0 = bestOuter.phi;
    else
        vmag0 = bestOuter.vmag;
        theta0 = bestOuter.theta;
        phi0 = bestOuter.phi;
    end

    % Use the same cost wrapper as the grid+fmin branch.
    % The bounds vmag_min/vmag_max only clamp inside the wrapper
    vmag_min = min(config.outerVmag);
    vmag_max = max(config.outerVmag);
    in_plane = numel(config.outerPhi) == 1 && abs(config.outerPhi) < 1e-6;

    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 500, 'Display', 'off');

    if in_plane
        obj = @(x) wrapperResonantCostForFminsearch([x(1), x(2), 0], ...
            config.vinfi, config.vinfo, config.body, config.jd0, config.T_p, ...
            config.n_val, config.m_val, config.apsis_flag, config.mu_sun, ...
            config.mu_planet, config.vmr_safety, vmag_min, vmag_max, ...
            config.res_flag, config.search_nu);
        x_opt = fminsearch(obj, [vmag0, theta0], options);
        vmag_opt  = max(vmag_min, min(vmag_max, x_opt(1)));
        theta_opt = x_opt(2);
        phi_opt   = 0;
    else
        obj = @(x) wrapperResonantCostForFminsearch(x, ...
            config.vinfi, config.vinfo, config.body, config.jd0, config.T_p, ...
            config.n_val, config.m_val, config.apsis_flag, config.mu_sun, ...
            config.mu_planet, config.vmr_safety, vmag_min, vmag_max, ...
            config.res_flag, config.search_nu);
        x_opt = fminsearch(obj, [vmag0, theta0, phi0], options);
        vmag_opt  = max(vmag_min, min(vmag_max, x_opt(1)));
        theta_opt = x_opt(2);
        phi_opt   = x_opt(3);
    end

    % Rebuild cartesian vinf_out and evaluate full breakdown
    vinf_out_refined = vmag_opt * [cos(phi_opt)*cos(theta_opt), cos(phi_opt)*sin(theta_opt), sin(phi_opt)];

    [dV_total, dV_GA1, dV_DSM, dV_GA2, ~, ~, ~, ~, ~, rp1, rp2] = evaluateTotalResonantCost(vinf_out_refined, config.vinfi, config.vinfo, config.body, config.jd0, config.T_p, config.n_val, config.m_val, config.apsis_flag, config.mu_sun, config.mu_planet, config.vmr_safety, config.res_flag, config.search_nu, true);

    if ~isfinite(dV_total)
        return;
    end

    % Refined nu via search_nu for the refined vinf_out
    [~, nu_best, revs_best, lp_best] = findOptimalDSMParameters( vinf_out_refined, config.body, config.jd0, config.T_p, config.n_val, config.m_val, config.apsis_flag, config.mu_sun, config.search_nu, true);

    outerSel = struct('theta', theta_opt, 'phi', phi_opt, 'vmag', vmag_opt, 'totalDv', dV_total, 'isFeasible', true, 'dV_GA1', dV_GA1, 'dV_DSM', dV_DSM, 'dV_GA2', dV_GA2, 'rp1', rp1, 'rp2', rp2, 'vinf_out', vinf_out_refined,'vinfX', vinf_out_refined(1), 'vinfY', vinf_out_refined(2), 'vinfZ', vinf_out_refined(3));
    innerSel = struct('nu', nu_best, 'totalDv', dV_total,'dV_DSM', dV_DSM, 'revs_before', revs_best, 'lp', lp_best, 'isFeasible', true);
end
