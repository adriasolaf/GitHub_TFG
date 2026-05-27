function [ key ] = aga_aux_duplicate_key ( ind, tol )
%aga_aux_duplicate_key Duplicate detection key
%   Quantizes one chromosome for near-duplicate detection.
%
% Inputs:
%   ind: AGA individual
%   tol: quantization tolerance per gene
%
% Outputs:
%   key: quantized duplicate detection key
%
% Example:
%   [ key ] = aga_aux_duplicate_key ( ind, tol );
%
% References:
%   [-]
%
%May 2026

    ind = ind(:)';
    scale = tol(:)';
    scale(scale == 0) = 1;
    key = round(ind ./ scale);
end
