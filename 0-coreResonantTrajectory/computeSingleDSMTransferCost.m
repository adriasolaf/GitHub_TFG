function [dV_DSM, r_m, v_m_minus, v_m_plus, vinf_in, va] = computeSingleDSMTransferCost(nu_DSM, revs_before, lp, N, apsis_flag, mu_sun, full_output_flag, orb_init, planets_state)
%   Constructs a single resonant orbit + DSM transfer for a specified
%   set of inputs. Given a DSM true anomaly, a revolution split, and a Lambert 
%   branch (lp), it validates the geometry via computeDSMStateAndTiming, solves a single 
%   Lambert arc to the arrival, and returns the resulting DSM Delta-V.
%
% Inputs:
%   nu_DSM: scanned angle for DSM true anomaly placement [rad]
%   revs_before: number of full Kepler revolutions before the DSM
%   lp: Lambert long-period flag (0 = short period, 1 = long period)
%   N: number of spacecraft revolutions in the resonant orbit
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: central body gravitational parameter [km^3/s^2]
%   full_output_flag: true = compute all output vectors,
%                     false = compute only dV_DSM
%
% Outputs:
%   dV_DSM: DSM Delta-V magnitude [km/s].
%   r_m: heliocentric position vector at the DSM point [km]
%   v_m_minus: spacecraft heliocentric velocity just before the DSM [km/s]
%   v_m_plus: spacecraft heliocentric velocity just after the DSM [km/s]
%   vinf_in: inbound hyperbolic excess velocity vector at arrival [km/s]
%   va: spacecraft heliocentric velocity at arrival [km/s]
%
% References:
%   [1] PhD thesis, David de la Torre Sangra
%   [2] Trajectory design of Solar Orbiter, J.M. Sanchez Perez
%   [3] The V-Infinity Leveraging Boundary Value Problem and Application in
%   Spacecraft Trajectory Design, Demyan V. Lantukh
%
% See also:
%   computeDSMStateAndTiming, findOptimalDSMParameters, Lambert, GetBodyICF
%
% Adria Sola Foixench
% April 2026

    % Initialize outputs
    dV_DSM = Inf; r_m = [NaN NaN NaN]; v_m_minus = [NaN NaN NaN]; v_m_plus = [NaN NaN NaN]; vinf_in = [NaN NaN NaN]; va = [NaN NaN NaN];

    % Validate geometry and get DSM state
    [flag_Res, r_m_cv, v_m_minus_cv, tof_lambert, ~] = computeDSMStateAndTiming(nu_DSM, revs_before, apsis_flag, mu_sun, orb_init, planets_state);

    if flag_Res ~= 0 
        return;
    end

    % Multi-revolution Lambert arc from DSM to arrival?
    mr_lambert = N - revs_before - 1;

    % Solve the single Lambert problem for the specified lp
    [v_m_plus_L, v_scf, flag_L, ~] = Lambert(r_m_cv, planets_state.r_pf, tof_lambert, mu_sun, [], mr_lambert, lp);

    % Check Lambert convergence
    if flag_L ~= 0
        return;
    end

    % Compute DSM Delta-V
    dV_DSM = norm(v_m_plus_L - v_m_minus_cv);
    r_m = r_m_cv;
    v_m_minus = v_m_minus_cv;

    % Full output
    if full_output_flag
        v_m_plus = v_m_plus_L;
        vinf_in = v_scf - planets_state.v_pf;
        va = v_scf;
    end

end
