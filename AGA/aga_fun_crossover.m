function [ child ] = aga_fun_crossover ( parentA, parentB, ~, ~, ranges )
%aga_fun_crossover Crossover resonant AGA chromosomes
%   Averages continuous genes and inherits the binary lp gene from one
%   randomly selected parent.
%
% Inputs:
%   parentA: first parent individual
%   parentB: second parent individual
%   ranges: chromosome bounds, one [min,max] row per gene
%
% Outputs:
%   child: descendant AGA individual
%
% Example:
%   [ child ] = aga_fun_crossover ( parentA, parentB, fitA, fitB, ranges );
%
% References:
%   [-]
%
%May 2026

    if nargin < 5 || isempty(ranges)
        ranges = aga_aux_default_ranges();
    end

    % Decode both parents so crossover receives bounded numeric vectors.
    a = aga_fun_decode(parentA, ranges);
    b = aga_fun_decode(parentB, ranges);

    % Average continuous parameters for a simple deterministic child.
    child = 0.5 * (a + b);

    % Keep lp binary by copying it from one parent.
    if rand() < 0.5
        child(5) = a(5);
    else
        child(5) = b(5);
    end

    child = aga_fun_decode(child, ranges);
end
