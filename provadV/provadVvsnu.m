body = 'Earth';
jd0 = -3281;
N = 1; 
M = 2; 
apsis_flag = 1; 
mu_sun = GetBodyProps('Sun');

vmag = 7.5;
theta_opt = deg2rad(225); 
phi_opt = deg2rad(-40);

vinf_guanyador = vmag * [cos(phi_opt)*cos(theta_opt), cos(phi_opt)*sin(theta_opt), sin(phi_opt)];

dV_vs_nu(body, jd0, vinf_guanyador, N, M, apsis_flag, mu_sun);