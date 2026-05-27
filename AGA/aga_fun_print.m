function aga_fun_print ( ind )
%aga_fun_print Print resonant AGA chromosome
%   Prints one AGA individual in a compact one-line format.
%
% Inputs:
%   ind: AGA individual [vinfo_km_s,theta_deg,phi_deg,nu_deg,lp]
%
% Outputs:
%   [-]
%
% Example:
%   aga_fun_print ( ind );
%
% References:
%   [-]
%
%May 2026

    fprintf('vinfo=%8.4f theta=%8.3f phi=%8.3f nu=%8.3f lp=%d ', ...
        ind(1), ind(2), ind(3), ind(4), round(ind(5)));
end
