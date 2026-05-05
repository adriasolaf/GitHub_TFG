function [r_arc1, r_arc2, r_m_out] = generateResonantTrajectoryPoints(r_p0, v_sc0, r_m, v_m_minus, v_m_plus, r_pf, va, revs_before_m, mr_lambert, mu_sun)
%   Calculates the sequence of position vectors for the complete
%   spacecraft trajectory, split into two arcs: from departure to the
%   DSM, and from the DSM to planetary arrival.
%
% Inputs:
%   r_p0: heliocentric position vector at departure [km]
%   v_sc0: spacecraft heliocentric velocity at departure [km/s]
%   r_m: heliocentric position vector at the DSM point [km]
%   v_m_minus: spacecraft velocity just before the DSM [km/s]
%   v_m_plus: spacecraft velocity just after the DSM [km/s]
%   r_pf: heliocentric position vector at arrival [km]
%   va: spacecraft heliocentric velocity at arrival [km/s]
%   revs_before_m: number of full revolutions before the DSM
%   mr_lambert: number of full revolutions in the Lambert arc
%   mu_sun: central body gravitational parameter [km^3/s^2]
%
% Outputs:
%   r_arc1: trajectory position points from departure to DSM [km] (Nx3)
%   r_arc2: trajectory position points from DSM to arrival [km] (Nx3)
%   r_m_out: exact position vector of the DSM point [km]
%
% References:
%   [-] n/a
%
% See also:
%   ICF2KEP_O, KEP2Arc, computeSingleDSMTransferCost
%
% Adria Sola Foixench
% April 2026

    % Arc 1: departure to DSM
    [sma1, ecc1, inc1, nu0_1, argp1, raan1] = ICF2KEP_O(r_p0, v_sc0, mu_sun);
    [~, ~, ~, nu_m_1, ~, ~] = ICF2KEP_O(r_m, v_m_minus, mu_sun);

    % Ensure forward propagation direction
    if nu_m_1 < nu0_1
        nu_m_1 = nu_m_1 + 2*pi;
    end
    nu_m_1_fullrev = nu_m_1 + (revs_before_m * 2*pi);

    % Generate arc 1 position points
    pts_arc1 = max(200, 200*(revs_before_m + 1));
    nu_vec_1 = linspace(nu0_1, nu_m_1_fullrev, pts_arc1);
    [r_arc1, ~] = KEP2Arc(sma1, ecc1, inc1, nu_vec_1, argp1, raan1, mu_sun);

    % Arc 2: DSM to arrival
    [sma2, ecc2, inc2, nu_m_2, argp2, raan2] = ICF2KEP_O(r_m, v_m_plus, mu_sun);
    [~, ~, ~, nuf_2, ~, ~] = ICF2KEP_O(r_pf, va, mu_sun);

    % Ensure forward propagation direction
    if nuf_2 < nu_m_2
        nuf_2 = nuf_2 + 2*pi;
    end
    nuf_2_fullrev = nuf_2 + (mr_lambert * 2*pi);

    % Generate arc 2 position points
    pts_arc2 = max(200, 200*(mr_lambert + 1));
    nu_vec_2 = linspace(nu_m_2, nuf_2_fullrev, pts_arc2);
    [r_arc2, ~] = KEP2Arc(sma2, ecc2, inc2, nu_vec_2, argp2, raan2, mu_sun);

    % Store DSM position
    r_m_out = r_m;

end
