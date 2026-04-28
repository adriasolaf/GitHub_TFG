function [ vinj ] = Vinjection ( vinf, rp, mu )
%VINJECTION Injection velocity
%   Returns the spacecraft velocity magnitude after an interplanetary orbit
%   injection manoeuvre
%
% Inputs:
%   vinf: spacecraft's planetocentric velocity at escape [km s^-1]
%   rp: initial radius where the injection manoeuvre is performed [km]
%   mu: standard gravitational parameter of the central body [km^3 s^-2]
%   
% Outputs:
%   vinj: spacecraft post-injection velocity magnitude [km s^-1]
%
% Example:
%   vinj = Vinjection ( 7, 300+6371, 4E5 ); % vinf 7km/s @ LEO 300km
%
% References:
%   [-]
%
%David de la Torre Sangra
%August 2014

% Escape velocity squared
% Vis-viva integral: vesc^2/2 - mu/rp = 0
vesc2 = 2.0 * mu ./ rp;

% Injection velocity
% Vis-viva integral: vinj^2/2 - mu/rp = vinf2/2 - 0; mu/rp = vesc^2/2
vinj = sqrt(vinf.^2 + vesc2);

end

