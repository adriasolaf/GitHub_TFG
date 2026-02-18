function [ dv, rp ] = GA_PGA2_Rp ( vinfi, vinfo, delta, mu  )
%GA_PGA2_Rp Powered gravity assist 2D
%   Computes a 2D powered gravity assist manoeuvre
%   Adjusts periapsis radius rp to match vinf_i/o and delta
%
% Inputs:
%   vinfi: spacecraft's planetocentric incoming velocity norm [km s^-1]
%   vinfo: spacecraft's planetocentric outgoing velocity norm [km s^-1]
%   delta: gravity assist turning angle [rad]
%   mu: standard gravitational parameter of the central body [km^3 s^-2]
%   
% Outputs:
%   dv: spacecraft's total velocity change norm [km s^-1]
%   rp: closest approach distance (periapsis distance) [km]
%
% Example:
%   [ dv, rp ] = GA_PGA2_Rp ( 3, 5, 0.4*pi, 4E5 );
%
% References:
%	[-]
%
%David de la Torre Sangra
%September 2018

% Parameters
ni = 1E2; % Maximum number of allowed iterations
tol = 1E-6; % Iteration tolerance

% Compute magnitudes of interest
a1 = mu / vinfi^2; % Semi-major axis of the arrival hyperbola
a2 = mu / vinfo^2; % Semi-major axis of the departure hyperbola

% Initial guess
% rp = a1; % Old version (sma of incoming hyperbola). Not very good.
% Mean point between the solution of each branch
rp1 = a1 * (1/sin(delta/2) - 1);
rp2 = a2 * (1/sin(delta/2) - 1);
rp = 0.5 * (rp1 + rp2);

% Iterate until convergence
for i=1:ni
    
    % Half-turning angle of incoming hyperbola
    delta1half = asin(a1 / (a1 + rp));
    
    % Half-turning angle of outgoing hyperbola
    delta2half = asin(a2 / (a2 + rp));
    
    % Turning angle of current iteration
    deltai = delta1half + delta2half;
    
    % Difference between current and target turning angles
    f = deltai - delta;
    
    % Derivative of incoming hyperbola's half-turning angle
    ddelta1half = - a1 / (sqrt(rp^2 + 2*a1*rp)*(a1 + rp));
    
    % Derivative of outgoing hyperbola's half-turning angle
    ddelta2half = - a2 / (sqrt(rp^2 + 2*a2*rp)*(a2 + rp));
    
    % Derivative of the difference between turning angles
    df = ddelta1half + ddelta2half;
    
    % Newton-Raphson iterative scheme on the periapsis radius
    rp_new = rp - f / df;
    
    % Safety checks
    if rp_new < 0 % Physically impossible
        rp = 0.5 * rp; % Halve the value of rp
        continue; % Keep iterating
    end
    
    % Check convergence
    if abs(rp_new - rp)/rp < tol % Convergence achieved
        break; % Stop iterating
    end
    
    % Update iterator value
    rp = rp_new;
    
end

% Compute deltaV (patching manoeuvre)
dv = abs(Vinjection(vinfo,rp,mu) - Vinjection(vinfi,rp,mu));

end

