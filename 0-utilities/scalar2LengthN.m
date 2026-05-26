function [vec, flag] = scalar2LengthN(x, n)
%   Expands a scalar input to a column vector of length n, or validates
%   that a vector input already has the required length.
%
% Inputs:
%   x: scalar or vector of length n
%   n: required output length (number of transfer legs)
%
% Outputs:
%   vec: column vector of length n (empty when flag ~= 0)
%   flag: 0 = success; 3 = invalid input shape (not scalar nor vector of length n)
%
% References:
%   [-] n/a
%
% See also:
%   MGA2_PGA2
%
% Adria Sola Foixench
% May 2026
    
    flag = 0;

    if isscalar(x)
        vec = repmat(x, n, 1);
    elseif isvector(x) && numel(x) == n
        vec = x(:);
    else
        vec = [];
        flag = 3;
    end
end
