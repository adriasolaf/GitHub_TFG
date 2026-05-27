function [ popOut, fitOut ] = aga_fun_unique ( popIn, fitIn, tol )
%aga_fun_unique Unique resonant AGA population
%   Removes near-duplicate chromosomes while keeping the first occurrence.
%
% Inputs:
%   popIn: input population cell array
%   fitIn: input fitness vector
%   tol: duplicate quantization tolerance per gene
%
% Outputs:
%   popOut: compact unique population cell array
%   fitOut: compact fitness vector
%
% Example:
%   [ popOut, fitOut ] = aga_fun_unique ( popIn, fitIn, tol );
%
% References:
%   [-]
%
%May 2026

    if isempty(popIn)
        popOut = popIn;
        fitOut = fitIn;
        return;
    end

    if nargin < 3 || isempty(tol)
        tol = [1e-6, 1e-5, 1e-5, 1e-5, 1];
    end

    seen = [];
    popOut = {};
    fitOut = [];
    for idx = 1:numel(popIn)
        if isempty(popIn{idx})
            continue;
        end

        % Quantize the chromosome so tiny floating-point differences do not
        % create separate individuals.
        key = aga_aux_duplicate_key(popIn{idx}, tol);

        % Append only the first copy of each quantized key.
        if isempty(seen) || ~any(all(abs(seen - key) < 1e-12, 2))
            seen(end + 1, :) = key; %#ok<AGROW>
            popOut{end + 1} = popIn{idx}; %#ok<AGROW>
            fitOut(end + 1, 1) = aga_aux_fitness_at(fitIn, idx); %#ok<AGROW>
        end
    end
end
