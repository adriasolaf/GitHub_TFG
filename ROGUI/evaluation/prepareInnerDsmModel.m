function [model, status] = prepareInnerDsmModel(config, vinf_out)
%prepareInnerDsmModel Build DSM search state for a fixed outgoing vinf.
%   This mirrors the setup section of findOptimalDSMParameters so GUI
%   sweeps can expose each nu/revolution/Lambert branch explicitly.
%
%   The output model contains the initial heliocentric conic, final planet
%   state after M planetary revolutions, and the compact structs expected by
%   computeSingleDSMTransferCost. It is shared by outer fixed-branch logic,
%   inner diagnostic maps, and trajectory reconstruction.

    status = struct('ok', false, 'message', '');
    model = struct();

    if nargin < 2 || numel(vinf_out) ~= 3 || any(~isfinite(vinf_out))
        status.message = 'Invalid outgoing v-infinity vector.';
        return;
    end

    % The resonant arrival epoch is defined by M revolutions of the
    % resonant body, not by an arbitrary GUI TOF for the DSM split.
    sec2days = 1 / 86400;
    model.tof_total_s = config.m_val * config.T_p;
    model.jd_f = config.jd0 + model.tof_total_s * sec2days;

    [model.r_p0, model.v_p0] = GetBodyICF(config.body, config.jd0, config.mu_sun, 0);
    [model.r_pf, model.v_pf] = GetBodyICF(config.body, model.jd_f, config.mu_sun, 0);

    model.planets_state.r_pf = model.r_pf;
    model.planets_state.v_pf = model.v_pf;
    model.planets_state.jd_f = model.jd_f;
    model.planets_state.tof_total_s = model.tof_total_s;

    model.vinf_out = vinf_out(:)';
    model.v_sc0 = model.v_p0 + model.vinf_out;

    % The pre-DSM arc must be an elliptic heliocentric orbit for the
    % resonance model used by computeDSMStateAndTiming.
    [sma, ecc, inc, nu0, argp, raan] = ICF2KEP_O(model.r_p0, model.v_sc0, config.mu_sun);
    if ecc >= 1.0 || sma <= 0
        status.message = 'Candidate produces a hyperbolic or invalid heliocentric initial orbit.';
        return;
    end

    n_motion = sqrt(config.mu_sun / sma^3);
    E0 = 2.0 * atan(sqrt((1.0 - ecc) / (1.0 + ecc)) * tan(nu0 / 2.0));
    M0 = E0 - ecc * sin(E0);

    model.orb_init.sma = sma;
    model.orb_init.ecc = ecc;
    model.orb_init.inc = inc;
    model.orb_init.argp = argp;
    model.orb_init.raan = raan;
    model.orb_init.n_motion = n_motion;
    model.orb_init.M0 = M0;
    model.orb_init.v_sc0 = model.v_sc0;

    status.ok = true;
    status.message = 'OK';
end
