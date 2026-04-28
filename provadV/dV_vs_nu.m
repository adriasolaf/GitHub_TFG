function dV_vs_nu(body, jd0, vinf_out, N, M, apsis_flag, mu_sun)
%   Plots dV_DSM vs nu search angle and Lambert arc eccentricity vs nu
%   for all Lambert branches, given a fixed VILM configuration.
%
% Inputs:
%   body: planet identifier string
%   jd0: departure epoch
%   vinf_out: outgoing v-infinity vector at departure [km/s]
%   N: number of spacecraft revolutions
%   M: number of body revolutions
%   apsis_flag: DSM reference apsis (1 = apoapsis, 0 = periapsis)
%   mu_sun: solar gravitational parameter [km^3/s^2]
%
% Adria Sola Foixench
% April 2026

    % Orbital period of the resonant body
    [sma_planet, ~, ~, ~, ~, ~] = GetBodyKEP_SSDG(body, jd0);
    T_p = 2*pi*sqrt(sma_planet^3 / mu_sun);

    % Time of flight and planet states
    sec2days = 1 / 86400;
    tof_total_s = M * T_p;
    jd_f = jd0 + (tof_total_s * sec2days);

    [r_p0, v_p0] = GetBodyICF(body, jd0, mu_sun, 0);
    [r_pf, v_pf] = GetBodyICF(body, jd_f, mu_sun, 0);

    planets_state.r_pf = r_pf;
    planets_state.v_pf = v_pf;
    planets_state.jd_f = jd_f;
    planets_state.tof_total_s = tof_total_s;

    % Spacecraft initial orbit
    v_sc0 = v_p0 + vinf_out;
    [sma, ecc_init, inc, nu0, argp, raan] = ICF2KEP_O(r_p0, v_sc0, mu_sun);

    if ecc_init >= 1.0
        return
    end

    n_motion = sqrt(mu_sun / (sma^3));
    E0 = 2.0 * atan(sqrt((1.0 - ecc_init) / (1.0 + ecc_init)) * tan(nu0 / 2.0));
    M0 = E0 - ecc_init * sin(E0);

    orb_init.sma = sma; orb_init.ecc = ecc_init; orb_init.inc = inc;
    orb_init.argp = argp; orb_init.raan = raan;
    orb_init.n_motion = n_motion; orb_init.M0 = M0;
    orb_init.v_sc0 = v_sc0;

    % Search grid
    revs_before = 0;
    nu_vec = linspace(0, pi - 1e-3, 1000);

    % Lambert branches
    mr_lambert_check = N - revs_before - 1;
    if mr_lambert_check > 0
        n_branches = 4;
        name_branches = {'lw=0, lp=0', 'lw=0, lp=1', 'lw=1, lp=0', 'lw=1, lp=1'};
        colors = {'b', 'g', 'r', 'k'};
    else
        n_branches = 2;
        name_branches = {'lw=0', 'lw=1'};
        colors = {'b', 'r'};
    end

    dV_mat = inf(n_branches, length(nu_vec));
    ecc_mat = inf(n_branches, length(nu_vec));

    % Evaluate all branches
    mr_lambert = N - revs_before - 1;
    
    idx_branch = 1;
    for lw = 0:1
        lp_max = (mr_lambert > 0) * 1;
        for lp = 0:lp_max
            for i = 1:length(nu_vec)
                [dV_DSM, r_m, ~, v_m_plus, ~, ~] = build_RessOrbDSM(nu_vec(i), revs_before, lw, lp, N, apsis_flag, mu_sun, true, orb_init, planets_state);
                dV_mat(idx_branch, i) = dV_DSM;
                if ~isinf(dV_DSM)
                    [~, ecc_lambert, ~, ~, ~, ~] = ICF2KEP_O(r_m, v_m_plus, mu_sun);
                    ecc_mat(idx_branch, i) = ecc_lambert;
                end
            end
            idx_branch = idx_branch + 1;
        end
    end

    % Departure anomaly in plot coordinates
    nu0_deg = rad2deg(nu0);
    nu_dep_deg = 180 - nu0_deg;
    if nu_dep_deg < 0
        nu_dep_deg = nu_dep_deg + 360;
    end

    % Plot
    figure('Name', sprintf('dV vs nu — %s N=%d M=%d', body, N, M), 'Color', 'w');

    % Upper subplot: Delta-V
    subplot(2, 1, 1);
    hold on;
    for k = 1:n_branches
        plot(rad2deg(nu_vec), dV_mat(k, :), 'Color', colors{k}, 'LineWidth', 1.5, 'DisplayName', name_branches{k});
    end
    xline(nu_dep_deg, 'm--', 'Sortida (\nu_0)', 'LineWidth', 1.5, 'FontSize', 10, ...
          'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    hold off;
    ylabel('\Delta V_{DSM} (km/s)', 'FontWeight', 'bold');
    title(sprintf('%s  N=%d M=%d  |  |v_\\infty|=%.2f km/s', body, N, M, norm(vinf_out)));
    grid on;
    legend('Location', 'best');
    finite_vals = dV_mat(~isinf(dV_mat));
    if ~isempty(finite_vals)
        ylim([max(0, min(finite_vals) - 0.5), min(finite_vals) + 10]);
    end

    % Lower subplot: eccentricity
    subplot(2, 1, 2);
    hold on;
    for k = 1:n_branches
        plot(rad2deg(nu_vec), ecc_mat(k, :), 'Color', colors{k}, 'LineWidth', 1.5);
    end
    yline(1.0, 'r--', 'Límit parabòlic (e=1)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    xline(nu_dep_deg, 'm--', 'Sortida (\nu_0)', 'LineWidth', 1.5, 'FontSize', 10, ...
          'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    hold off;
    xlabel('Angle de cerca \nu (deg)', 'FontWeight', 'bold');
    ylabel('Excentricitat arc Lambert', 'FontWeight', 'bold');
    grid on;
    ylim([0, 1.5]);

end
