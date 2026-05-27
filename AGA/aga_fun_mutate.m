function [ mutant ] = aga_fun_mutate ( ind, ~, ranges, mutationStep )
%aga_fun_mutate Mutate resonant AGA chromosome
%   Perturbs continuous genes with Gaussian noise and occasionally flips
%   the binary Lambert long-period flag.
%
% Inputs:
%   ind: AGA individual [vinfo_km_s,theta_deg,phi_deg,nu_deg,lp]
%   ranges: chromosome bounds, one [min,max] row per gene
%   mutationStep: standard deviations for [vinfo,theta,phi,nu]
%
% Outputs:
%   mutant: mutated AGA individual
%
% Example:
%   [ mutant ] = aga_fun_mutate ( ind, fit, ranges, mutationStep );
%
% References:
%   [-]
%
%May 2026

    if nargin < 3 || isempty(ranges)
        ranges = aga_aux_default_ranges();
    end
    if nargin < 4 || isempty(mutationStep)
        mutationStep = [1.0, 20.0, 20.0, 20.0];
    end
    if isscalar(mutationStep)
        mutationStep = repmat(mutationStep, 1, 4);
    end

    % Start from a valid chromosome so mutation never propagates bad shape.
    mutant = aga_fun_decode(ind, ranges);

    % Perturb [vinfo, theta, phi, nu] with caller-controlled step sizes.
    mutant(1:4) = mutant(1:4) + mutationStep(1:4) .* randn(1, 4);

    % Flip lp occasionally to preserve exploration of both Lambert branches.
    if rand() < 0.15
        mutant(5) = 1 - mutant(5);
    end

    % Clamp final genes after mutation.
    mutant = aga_fun_decode(mutant, ranges);
end
