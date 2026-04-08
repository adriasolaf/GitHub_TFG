function current_valid = checkLambertValidity(r_m_out, r_pf, tof_lambert_out, mu_sun, mr_lambert_try)
%   Evaluates if a simple/multi-revolution Lambert boundary value problem
%   is geometrically possible and if the available TOF
%   exceeds the minimum TOF required for the specified revolutions.
%
% Inputs:
%   r_m_out: heliocentric position vector at the DSM point [km]
%   r_pf: heliocentric position vector at the final arrival point [km]
%   tof_lambert_out: time of flight available for the Lambert transfer arc [s]
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   mr_lambert_try: number of full spacecraft revolutions for the Lambert solver
%
% Outputs:
%   current_valid: boolean flag indicating if the Lambert problem is valid
%
% References:
%   [1] PhD_DavidDLTS_astro (David de la Torre)
%
% See also:
%   Lambert_Izzo_2015_X_Tmin
%
% Adrià Solà Foixench
% April 2026

    current_valid = false;

    % Check collinearity for non degenerate geometry
    cross_rm_rpf = cross(r_m_out, r_pf);
    sin_theta = norm(cross_rm_rpf) / (norm(r_m_out) * norm(r_pf));

    if sin_theta < 1e-3
        return; % Too close to collinear
    end

    % Normalize Lambert inputs
    [nr1, nr2, ntof, nmu, ~, ~, ~] = NorMag_Lambert(r_m_out, r_pf, tof_lambert_out, mu_sun);
    r1n = norm(nr1); % Normalized radius 1
    r2n = norm(nr2); % Normalized radius 2
    c = norm(nr2 - nr1); % Chord length
    s = 0.5*(r1n + r2n + c); % Semi-perimeter

    % Normalized time of flight
    t = sqrt(2*nmu / (s^3)) * ntof;

    % Compute minimum TOF for multi-rev Lambert
    [~, Tmin_norm_0, ~, ~] = Lambert_Izzo_2015_X_Tmin(nr1, nr2, nmu, 0, mr_lambert_try);
    [~, Tmin_norm_1, ~, ~] = Lambert_Izzo_2015_X_Tmin(nr1, nr2, nmu, 1, mr_lambert_try);
    Tmin_0 = Tmin_norm_0 * sqrt(2*nmu / (s^3)); % Minimum TOF way 0
    Tmin_1 = Tmin_norm_1 * sqrt(2*nmu / (s^3)); % Minimum TOF way 1

    % Valid only if available TOF exceeds minimum for at least one way
    if t > Tmin_0 || t > Tmin_1
        current_valid = true;
    end

end