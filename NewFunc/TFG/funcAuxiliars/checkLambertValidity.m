function current_valid = checkLambertValidity(r_m_out, r_pf, tof_lambert_out, mu_sun, mr_lambert_try)
%   Evaluates if a multi-revolution Lambert problem is
%   geometrically possible and if the available TOF exceeds the minimum
%   TOF required for the specified number of revolutions.
%
% Inputs:
%   r_m_out: heliocentric position vector at the DSM point [km]
%   r_pf: heliocentric position vector at the arrival point [km]
%   tof_lambert_out: time of flight available for the Lambert arc [s]
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   mr_lambert_try: number of full revolutions for the Lambert solver
%
% Outputs:
%   current_valid: boolean, true if the Lambert problem is feasible
%
% References:
%   [1] PhD thesis, David de la Torre Sangra
%
% See also:
%   Lambert_Izzo_2015_X_Tmin, build_RessOrbDSM
%
% Adria Sola Foixench
% April 2026

    current_valid = false;

    % Degenerate geometry if r_m and r_pf are nearly parallel
    cross_rm_rpf = cross(r_m_out, r_pf);
    sin_theta = norm(cross_rm_rpf) / (norm(r_m_out) * norm(r_pf));

    if sin_theta < 1e-3
        return; % Too close to collinear for a posible Lambert solution
    end

    % Normalize Lambert inputs for minimum TOF computation
    [nr1, nr2, ntof, nmu, ~, ~, ~] = NorMag_Lambert(r_m_out, r_pf, tof_lambert_out, mu_sun);
    r1n = norm(nr1);
    r2n = norm(nr2);
    c = norm(nr2 - nr1); % Chord length
    s = 0.5*(r1n + r2n + c); % Semi-perimeter

    % Normalized time of flight
    t = sqrt(2*nmu / (s^3)) * ntof;

    % Compute minimum TOF for both short-way and long-way transfers
    [~, Tmin_norm_0, ~, ~] = Lambert_Izzo_2015_X_Tmin(nr1, nr2, nmu, 0, mr_lambert_try);
    [~, Tmin_norm_1, ~, ~] = Lambert_Izzo_2015_X_Tmin(nr1, nr2, nmu, 1, mr_lambert_try);
    Tmin_0 = Tmin_norm_0 * sqrt(2*nmu / (s^3));
    Tmin_1 = Tmin_norm_1 * sqrt(2*nmu / (s^3));

    % Valid if available TOF exceeds minimum for at least one way
    if t > Tmin_0 || t > Tmin_1
        current_valid = true;
    end

end
