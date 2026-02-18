function [ dv, extra ] = GA_PGA2_Vinfo ( vinfi, vinfo, delta, rp, mu )
%GA_PGA2_Vinfo Powered gravity assist 2D
%   Computes a 2D powered gravity assist manoeuvre
%   Matches hyperbolas at at rp keeping delta, and then fixes vinfo
%
% Inputs:
%   vinfi: spacecraft's planetocentric incoming velocity norm [km s^-1]
%   vinfo: spacecraft's planetocentric outgoing velocity norm [km s^-1]
%   delta: gravity assist turning angle [rad]
%   mu: standard gravitational parameter of the central body [km^3 s^-2]
%
% Outputs:
%   dv: spacecraft's total velocity change norm [km s^-1]
%   extra: extra output struct
%
% Example:
%   [ dv, extra ] = GA_PGA2_Vinfo ( 3, 5, 0.4*pi, 4E3, 4E5 );
%
% References:
%	[-]
%
%David de la Torre Sangra
%September 2018

% Compute magnitudes of interest
a1 = mu / vinfi^2; % Semi-major axis of the arrival hyperbola

% Half-turning angle of incoming hyperbola
delta1half = asin(a1 / (a1 + rp));

% Half-turning angle of outgoing hyperbola
delta2half = delta - delta1half;

% Semi-major axis of the departure hyperbola
a2 = rp / (1/sin(delta2half) - 1);

% Vinf (virtual) of departure hyperbola
vinfov = sqrt(mu / a2);

% Compute deltaV (patching manoeuvre)
dv_rp = abs(Vinjection(vinfov,rp,mu) - Vinjection(vinfi,rp,mu));

% Fix vinfo modulus
dv_inf = abs(vinfo - vinfov);

% Total dv modulus
dv = dv_rp + dv_inf;

% Assemble extra inputs
extra.a1 = a1;
extra.a2 = a2;
extra.dv_rp = dv_rp;
extra.dv_inf = dv_inf;
extra.delta1half = delta1half;
extra.delta2half = delta2half;
extra.vinfov = vinfov;

end

