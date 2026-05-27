function count = costEvalCounter(action)
%costEvalCounter Persistent counter for resonant cost evaluations.
%   action can be 'inc', 'get', or 'reset'. The function is intentionally
%   lightweight so callers can ignore its output.

    persistent nEval
    if isempty(nEval)
        nEval = 0;
    end

    if nargin < 1 || isempty(action)
        action = 'get';
    end

    switch lower(action)
        case 'inc'
            nEval = nEval + 1;
        case 'reset'
            nEval = 0;
        case 'get'
            % No state change.
        otherwise
            error('costEvalCounter:InvalidAction', 'Unknown action: %s', action);
    end

    count = nEval;
end
