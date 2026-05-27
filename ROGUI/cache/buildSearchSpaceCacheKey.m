function key = buildSearchSpaceCacheKey(config, item)
%buildSearchSpaceCacheKey Build ROGUI cache key
%   Creates a deterministic string key from the mission configuration and
%   optional selected outer or inner search row.
%
% Inputs:
%   config: normalized ROGUI mission/search configuration structure
%   item: optional table row or structure with candidate variables
%
% Outputs:
%   key: deterministic cache key string
%
% Example:
%   [ key ] = buildSearchSpaceCacheKey ( config, selectedRow );
%
% References:
%   [-]
%
%May 2026

    mission = strjoin(config.planets, '-');
    base = sprintf('%s|jd%.6f|tof%s|N%s|M%s|leg%d', mission, config.jd2k0, ...
        encodeNumeric(config.tofs), encodeNumeric(config.N), encodeNumeric(config.M), config.resonantIndex);

    if nargin < 2 || isempty(item)
        key = base;
        return;
    end

    if istable(item)
        item = table2struct(item(1, :));
    end

    % Use rounded text rather than raw binary floating point. This is stable
    % enough for GUI cache reuse while still distinguishing neighboring
    % grid samples.
    parts = {base};
    names = {'vmag', 'theta', 'phi', 'nu', 'revs_before', 'lp', 'vinfX', 'vinfY', 'vinfZ'};
    for idx = 1:numel(names)
        if isfield(item, names{idx})
            parts{end + 1} = sprintf('%s=%.10g', names{idx}, item.(names{idx})); %#ok<AGROW>
        end
    end
    key = strjoin(parts, '|');
end

function txt = encodeNumeric(values)
    txt = sprintf('_%.10g', values(:));
end
