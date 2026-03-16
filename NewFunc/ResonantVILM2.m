function [ dV, vinf_in, v_sc0, va ] = ResonantVILM2( body, jd0, vinf_out, N, M, apsis_flag, mu_sun )

sec2days = 1 / 86400;


T_p = 365.25636 * 24 * 3600;
[ r_p0, v_p0 ] = GetBodyICF(body, jd0, mu_sun, 0);


tof_total_s = M * T_p; 
jd_f = jd0 + (tof_total_s * sec2days);

[ r_pf, v_pf ] = GetBodyICF(body, jd_f, mu_sun, 0);


v_sc0 = v_p0 + vinf_out;

[sma, ecc, inc, nu0, argp, raan ] = ICF2KEP_O( r_p0, v_sc0, mu_sun );

if ecc >= 1.0 % Perque no grafiqui orbites no objectiu 
    dV = NaN; 
    vinf_in = [NaN NaN NaN]; 
    va = [NaN NaN NaN]; 
    return;
end

n = sqrt(mu_sun / (sma^3));

E0 = 2.0 * atan(sqrt((1.0 - ecc) / (1.0 + ecc)) * tan(nu0 / 2.0));
M0 = E0 - ecc * sin(E0);

canv_angle = 0;        
step1 = 0.05;     
step2 = 0.00025;  
max_iter = 150;        
valid = false;
is_step1 = true;

for iter = 1:max_iter

    offset = 1e-3;
    if apsis_flag == 1
        nu_m = pi + offset - canv_angle; 
    else
        nu_m = 0 + offset - canv_angle; 
    end
    
    nu_m = mod(nu_m, 2*pi); % Limitar la volta dins de 2pi
    
    Em = 2.0 * atan(sqrt((1.0 - ecc) / (1.0 + ecc)) * tan(nu_m / 2.0));
    Mm = Em - ecc * sin(Em);
    
    delta_M = Mm - M0;

    if delta_M <= 0
        delta_M = delta_M + 2*pi; 
    end
    
    revs_before_m = floor(N/2);  % Agafa el numero enter de la divisio
    delta_M = delta_M + (revs_before_m * 2 * pi);
    mr_lambert = N - revs_before_m - 1; 
    
    tof_m_s = delta_M / n;   % Converitr angle a temps
    jd_m = jd0 + (tof_m_s * sec2days);
    
    if jd_m >= jd_f
        if is_step1
            canv_angle = canv_angle + step1;
        else
            canv_angle = canv_angle + step2;
        end
        continue;
    end
    
    [ r_m, v_m_minus ] = KEP2ICF_O( sma, ecc, inc, nu_m, argp, raan, mu_sun );
    
    
    tof_lambert_s = tof_total_s - tof_m_s;
    
   
    actual_valid = false;
    if mr_lambert > 0
        theta = DeltaNu3(r_m, r_pf, 1); 
        if theta > pi
            lw = 1; 
        else 
            lw = 0; 
        end 
        
        % Comprovacio extreta de Tesis David de la Torre
        [ nr1, nr2, ntof, nmu, ~, ~, ~ ] = NorMag_Lambert( r_m, r_pf, tof_lambert_s, mu_sun );
        r1n = norm(nr1); 
        r2n = norm(nr2);
        c = norm(nr2 - nr1); 
        s = 0.5 * (r1n + r2n + c);
        t = sqrt(2 * nmu / (s*s*s)) * ntof;
        
        
        [ ~, Tmin_norm, ~, ~ ] = Lambert_Izzo_2015_X_Tmin( nr1, nr2, nmu, lw, mr_lambert );
        Tmin = Tmin_norm * sqrt(2 * nmu / (s*s*s));
        
        if Tmin < t
            actual_valid = true;
        end
    else
        actual_valid = true;
    end
    
    if actual_valid
        if is_step1 && canv_angle > 0
            canv_angle = canv_angle - step1;
            is_step1 = false;
        else
            valid = true; 
            break;
        end
    else
        if is_step1
            canv_angle = canv_angle + step1;
        else
            canv_angle = canv_angle + step2;
        end
    end
end

if ~valid
    dV = NaN; 
    vinf_in = [NaN NaN NaN];
    va = [NaN NaN NaN];
    return;
end


dV = NaN; 
vinf_in = [NaN NaN NaN];
va = [NaN, NaN, NaN];
best_dV = inf;

lp_max = 0;

if mr_lambert > 0 
    lp_max = 1; 
end

for lp = 0:lp_max
    for lw = 0:1
        [ v_m_d, v_scf_a, flag, ~ ] = Lambert( r_m, r_pf, tof_lambert_s, mu_sun, lw, mr_lambert, lp );
    
        if flag == 0

            dV_tmp = norm(v_m_d - v_m_minus);
 
            if dV_tmp < best_dV
                best_dV = dV_tmp;
                dV = dV_tmp;
                vinf_in = v_scf_a - v_pf;
                va = v_scf_a;
            end
        end
    end
end

end