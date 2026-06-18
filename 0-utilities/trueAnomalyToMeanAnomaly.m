function M = trueAnomalyToMeanAnomaly(nu, ecc)
% trueAnomalyToMeanAnomaly Convert elliptic true anomaly to mean anomaly.
%
% Inputs:
%   nu: true anomaly [rad]
%   ecc: eccentricity, elliptic orbit only
%
% Outputs:
%   M: mean anomaly wrapped to [0, 2*pi) [rad]

    E = atan2(sqrt(1.0 - ecc^2) .* sin(nu), ecc + cos(nu));
    M = mod(E - ecc .* sin(E), 2*pi);
end
